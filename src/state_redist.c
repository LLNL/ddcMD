#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "integer.h"
#include "treecomm.h"
#include "balance.h"
#include "redistribute.h"

//#include "state.h"
//#include "species.h"
//#include "group.h"
#include "state_redist.h"

#include "ptiming.h"
#include "heap.h"


typedef struct {
  double xx[3],vv[3];
  gid_type label;
  int atomtype,species_idx,group_idx;
} state_aos;

integer redistribute_longcut(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			     integer nloc,integer nlocmax,state_aos *data,
			     double c0[3],double c1[3],cmp_fun cmp_coord[3]) {
  const int np = (pid1-pid0)/pstride;
  const int g_pid = (pid-pid0)/pstride;
  integer pcut;
  int idir;
  double ccut0[3],ccut1[3];

  if(np <= 1) return nloc;

  /* Find longest box dimension and cut along that */ {
    int i;
    double xcut;

    idir = 0;
    for(i = 1; i<3; i++) {
      if(c1[i]-c0[0] > c1[idir]-c0[idir])
	idir = i;
    }

    profile(RECBIS_BALANCE2,START);
    /* Cut along idir direction, and extract cutting pivot value */ {
      state_aos pivot;
      nloc = balance2(pid0,pid1,pid,pstride,comm,
		      nloc,nlocmax,sizeof(*data),data,cmp_coord[idir],&pivot);
      xcut = pivot.xx[idir];
    }
    profile(RECBIS_BALANCE2,END);

    /* Determine cutting plane */ {
      double x = xcut,dx = fabs(x),y = x+dx,y1;
      if(dx == 0.0) y = dx = 1.0;
      do {
	y1 = y;
	dx *= 0.5;
	y = x + dx;
      } while(y > x);
      for(i = 0; i<3; i++) {
	ccut0[i] = c0[i];
	ccut1[i] = c1[i];
      }
      ccut0[idir] = y1;
      ccut1[idir] = y1;
    }
  }
  
  pcut = (np+1)/2;
  if(g_pid < pcut) {
    nloc = redistribute_longcut(pid0,pid0+pcut*pstride,pid,pstride,comm,
				nloc,nlocmax,data,c0,ccut1,cmp_coord);
  } else {
    nloc = redistribute_longcut(pid0+pcut*pstride,pid1,pid,pstride,comm,
				nloc,nlocmax,data,ccut0,c1,cmp_coord);
  }

  return nloc;
}




static int cmp_coord(const void *ain,const void *bin,int icoord) {
 const state_aos
   *a = (const state_aos *) ain,
   *b = (const state_aos *) bin;
 double
   ax = a->xx[icoord],
  bx = b->xx[icoord];
 if(ax < bx) return -1;
 else if(ax > bx) return 1;
 else if(ax == bx) return 0;

 /* We can end up here e.g. if ax or bx is NAN. */
 assert(1 == 0);

 /* Nans compare false to anything, so... */ {
   union { double x; unsigned long long int i; } au,bu;
   assert(sizeof(unsigned long long int) >= sizeof(double));
   au.i = 0;
   bu.i = 0;
   au.x = ax;
   bu.x = bx;
   if(au.i < bu.i) return -1;
   else if(au.i > bu.i) return 1;
   else return 0;
 }
}

static int cmp_x(const void *ain,const void *bin) { return cmp_coord(ain,bin,0); }
static int cmp_y(const void *ain,const void *bin) { return cmp_coord(ain,bin,1); }
static int cmp_z(const void *ain,const void *bin) { return cmp_coord(ain,bin,2); }

static void accum_max_integer(void *x,const void *dx) {
  integer *ip = (integer *) x;
  const integer *dip = (const integer *) dx;
  if(dip[0] > ip[0]) ip[0] = dip[0];
}

integer state_redist(STATE *state_p,MPI_Comm comm,integer nitems) {
  const int debug = 0;
  integer n = nitems,bufsz;

  state_aos *atoms;

  int pid,np;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&pid);

  /* Determine needed buffer size... */ if(0) {
    integer nmax = nitems;
    tree_allreduce(0,np,pid,1,comm,
		   1,sizeof(integer),&nmax,accum_max_integer);
    bufsz = 3*nmax + 10000;

    if(debug > 0) {
      if(pid == 0)
	printf("{%s:%d} pid=%03d  in %s(), atom buffer allocated, nitems=%d, nmax=%d, bufsz=%d, np=%d\n",
	       __FILE__,__LINE__,pid,__func__,
	       (int) nitems,(int) nmax,(int) bufsz,np);
    }
  }


  //atoms = (state_aos *) malloc(sizeof(state_aos) * bufsz);
  unsigned int atomsBlock;
  atoms = (state_aos *) heapGet(&atomsBlock);
  assert(atoms != NULL && "Heap error for atoms" != NULL);
  bufsz = heapAvailable()/sizeof(*atoms);

  if(debug > 0) {
    tree_barrier(0,np,pid,1,comm);
    if(pid == 0)
      printf("{%s:%d} pid=%03d  in %s(), atom buffer allocated, n=%d, bufsz=%d, np=%d -- state_p->nlocal=%d\n",
	     __FILE__,__LINE__,pid,__func__,
	     (int) n,(int) bufsz,np,(int) (state_p->nlocal));
  }

  assert(bufsz >= n);
  for(integer i = 0; i<n; i++) {
    atoms[i].xx[0] = state_p->rx[i];
    atoms[i].xx[1] = state_p->ry[i];
    atoms[i].xx[2] = state_p->rz[i];
    atoms[i].vv[0] = state_p->vx[i];
    atoms[i].vv[1] = state_p->vy[i];
    atoms[i].vv[2] = state_p->vz[i];

    atoms[i].atomtype = state_p->atomtype[i];
    atoms[i].label = state_p->label[i];
    atoms[i].species_idx = state_p->species[i]->index;
    atoms[i].group_idx = state_p->group[i]->index;
  }
  
  if(debug > 0) {
    tree_barrier(0,np,pid,1,comm);
    if(pid == 0)
      printf("{%s:%d} pid=%03d  in %s(), data copied to atom buffer, calling redistribute()\n",
	     __FILE__,__LINE__,pid,__func__);
  }

  profile(RECBIS_REDIST,START);
  if(0) {
    /* Redistribute atoms to equalize... */   
    cmp_fun cmp_list[3] = { cmp_x , cmp_y , cmp_z };
    
    n = redistribute(0,np,pid,1,comm,
		     n,bufsz,sizeof(state_aos),atoms,3,0,cmp_list);
  } else {
    /* Redistribution which always cuts along longest dimension, and thus
       keeps domains cluser to cubic (spherical) */
    cmp_fun cmp_list[3] = { cmp_x , cmp_y , cmp_z };
    double
      c0[3] = { INFINITY, INFINITY, INFINITY},
      c1[3] = {-INFINITY,-INFINITY,-INFINITY},
      c0glob[3],c1glob[3];
    int i;
    for(i = 0; i<n; i++) {
      int j;
      for(j = 0; j<3; j++) {
	if(atoms[i].xx[j] < c0[j]) c0[j] = atoms[i].xx[j];
	if(atoms[i].xx[j] > c1[j]) c1[j] = atoms[i].xx[j];
      }
    }
    if(n == 0) {
      int j;
      for(j = 0; j<3; j++) {
	c0[j] =  INFINITY;
	c1[j] = -INFINITY;
      }
    }
    MPI_Allreduce(c0,c0glob,3,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(c1,c1glob,3,MPI_DOUBLE,MPI_MAX,comm);
    
    n = redistribute_longcut(0,np,pid,1,comm,
			     n,bufsz,atoms,
			     c0glob,c1glob,cmp_list);
  }
  profile(RECBIS_REDIST,END);

  if(debug > 0) {
    tree_barrier(0,np,pid,1,comm);
    if(pid == 0)
      printf("{%s:%d} pid=%03d  in %s(), new n = %d, calling resize()\n",
	     __FILE__,__LINE__,pid,__func__,(int) n);
  }

  resize(/* new size = */ n,/* mode = */ 2,state_p);
  state_p->nlocal = n;

  if(debug > 0) {
    tree_barrier(0,np,pid,1,comm);
    if(pid == 0)
      printf("{%s:%d} pid=%03d  in %s(), resize returned. Copying atom data back to state array.\n",
	     __FILE__,__LINE__,pid,__func__);
  }

  /* Copy redistributed atoms to simulation/state array */ {
    const double nth = 1.0/n;
    double xxavg[3] = { 0.0 , 0.0 , 0.0 };

    for(integer i = 0; i<n; i++) {
      int j;
      for(j = 0; j<3; j++)
	xxavg[j] += atoms[i].xx[j] * nth;

      state_p->rx[i] = atoms[i].xx[0];
      state_p->ry[i] = atoms[i].xx[1];
      state_p->rz[i] = atoms[i].xx[2];
      state_p->vx[i] = atoms[i].vv[0];
      state_p->vy[i] = atoms[i].vv[1];
      state_p->vz[i] = atoms[i].vv[2];
      
      state_p->atomtype[i] = atoms[i].atomtype;
      state_p->label[i] = atoms[i].label;
      state_p->species[i] = species_by_index(NULL,atoms[i].species_idx);
      state_p->group[i] = group_by_index(NULL,atoms[i].group_idx);
    }
    //free(atoms);
    heapEndBlock(atomsBlock,n * sizeof(*atoms));
    heapFree(atomsBlock);

    /* Communicate domain centers resulting from redistribution */ {
      double (*domcen)[3] = (double (*)[3]) malloc(sizeof(double [3]) * np);
      /* This may actually be done by Allgrather, but I'm too
	  lazy to perform semantics research/tests now. */ {
	int j;
	for(j = 0; j<3; j++) domcen[pid][j] = xxavg[j];
      }
      MPI_Allgather(xxavg,3,MPI_DOUBLE,domcen,3,MPI_DOUBLE,comm);

      state_p->domcen_active = 1;
      state_p->domcen_np = np;
      state_p->domcen_pid = pid;
      state_p->domcen_vec = domcen;

    }

  }

  if(debug > 0) {
    tree_barrier(0,np,pid,1,comm);
    if(pid == 0)
      printf("{%s:%d} pid=%03d  in %s(), copying complete. Returning.\n",
	     __FILE__,__LINE__,pid,__func__);
  }

  return n; 
}
