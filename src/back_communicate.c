#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <mpi.h>

#include "ddc.h"

#include "genptr.h"
#include "integer.h"
#include "treecomm.h"

#include "equalize.h"
#include "exchange.h"
#include "balance.h"

#include "rectimer.h"
#include "functimer.h"


#include "balance_hist.h"
#include "back_communicate.h"
#include "alltoall_sparse.h"
#include "simulate.h"
#include "bisectionLoadBalance.h"

#include "utilities.h"
#include "ptiming.h"
#include "heap.h"

#include "recnlist.h"

typedef struct {
  int origin,destination;
  int index_at_origin;
} back_comm_data;

integer get_destination(const void *ptr) {
  const back_comm_data *bcdp = (const back_comm_data *) ptr;
  return (integer) bcdp->origin;
}

#include "bisection_data.h"
/*
typedef struct {
  float coord[3];
  int origin,index_at_origin;
} bisection_data;
*/
int cmp_coord(const void *ap,const void *bp,int idx) {
  double
    a = ((bisection_data *) ap)->coord[idx],
    b = ((bisection_data *) bp)->coord[idx];
  if(a <  b) return -1;
  if(a >  b) return  1;
  if(a == b) return  0;
  /* Take care of nan's */ {
    return memcmp(&(((bisection_data *) ap)->coord[idx]),
		  &(((bisection_data *) bp)->coord[idx]),
		  sizeof( *(((bisection_data *) ap)->coord) ));
  }
}
int cmp_coord_x(const void *ap,const void *bp) { return cmp_coord(ap,bp,0); }
int cmp_coord_y(const void *ap,const void *bp) { return cmp_coord(ap,bp,1); }
int cmp_coord_z(const void *ap,const void *bp) { return cmp_coord(ap,bp,2); }

void back_communicate(integer pid,integer np,MPI_Comm comm,
		      integer nloc,integer nlocmax,bisection_data *bd_p,
		      PARTICLE *particle_array) {
  STARTTIMER;
  integer i;
  back_comm_data *bcd_p = (back_comm_data *) bd_p;
  const integer nlocmax_bcd = (nlocmax*sizeof(bisection_data)) / sizeof(back_comm_data);
  int nself = 0;

  //if(pid == 0) printf("Entering back_communicate...\n");
  assert(sizeof(bisection_data) >= sizeof(back_comm_data));

  for(i = 0; i<nloc; i++) {
    bisection_data tmp = bd_p[i];
    integer pos;
    if(tmp.origin == pid) {
      pos = nself++;
      if(i > pos)
	bcd_p[i] = bcd_p[pos];
    } else
      pos = i;
    bcd_p[pos].origin = tmp.origin;
    bcd_p[pos].destination = pid;
    bcd_p[pos].index_at_origin = tmp.index_at_origin;
  }

  profile(RECBIS_A2ASPARSE,START);
  nloc = nself + alltoall_sparse(0,np,pid,1,comm,
				 nloc-nself,nlocmax_bcd-nself,sizeof(back_comm_data),
				 bcd_p+nself,&get_destination);
  profile(RECBIS_A2ASPARSE,END);
  /* Debug return data, make sure each receive location is visited exactly once. */ if(0) {
    int *tag = (int *) malloc(sizeof(int) * nloc);
    integer i;
    for(i = 0; i<nloc; i++)
      tag[i] = 0;
    for(i = 0; i<nloc; i++) {
      const integer j = bcd_p[i].index_at_origin;
      if(bcd_p[i].origin != pid) {
	printf("(pid=%03d) Error at %s:%d in %s(): pid(%d) != origin(%d)...\n",
	       (int) pid,__FILE__,__LINE__,__func__,(int) pid,(int) bcd_p[i].origin);
	exit(1);
      }
      if(j <0 || j>= nloc) {
	printf("(pid=%03d) Error at %s:%d in %s(): "
	       "origin index j=%d out of bounds [0,%d)...\n",
	       (int) pid,__FILE__,__LINE__,__func__,(int) j,(int) nloc);
	exit(1);	
      }
      tag[j]++;
    }
    for(i = 0; i<nloc; i++) {
      if(tag[i] != 1) {
	printf("(pid=%03d) Error at %s:%d in %s(): "
	       "Counting error, tag[%d] = %d\n",
	       (int) pid, __FILE__,__LINE__,__func__,(int) i,tag[i]);
	exit(1);	
      }
    }
    free(tag);
  }

  for(i = 0; i<nloc; i++) {
    particle_array[bcd_p[i].index_at_origin].domain_id = bcd_p[i].destination;
  }
  STOPTIMER;
  //if(pid == 0) printf("leaving back_communicate\n");
}

#if 0
/* This routine actually may create a comm world sized array :-( */
integer treesize(integer n) {
  static integer sz = 0,*list = NULL;

  if(sz < n+1) {
    integer nsz = sz*2,i;
    if(nsz < n+1) nsz = n+1;
    list = (integer *) realloc(list,sizeof(integer) * nsz);
    assert(list != NULL);
    for(i = sz; i<nsz; i++)
      list[i] = (i<2) ? i : -1;
    sz = nsz;
  }

  if(list[n] < 0) {
    integer n2 = (n+1)/2;
    list[n] = 1 + treesize(n2) + treesize(n-n2);
  }
  return list[n];
}
#else
/* However, this one runs on O(log n) and uses constant storage :-) */
integer treesize(integer n) {
  integer s = 0,a = 1,b = 0;
  if(n == 0) return 0;
  while(n > 1) {
    s += a + b;
    if(n & 1)
      b = a + 2*b;
    else
      a = 2*a + b;
    n = (n+1) >> 1;
  }
  return s + a - b;
}
#endif

integer treeself(integer p0,integer p1,integer p) {
  integer pcut = p0 + (p1-p0+1)/2;
  if(p1 - p0 <= 1) return 0;
  else if(p < pcut) return 1 + treeself(p0,pcut,p);
  else return 1 + treesize(pcut) + treeself(0,p1-pcut,p-pcut);
}

void min_double(void *min_tot,const void *test_val) {
  double dx = *(double *) test_val;
  double *xp = (double *) min_tot;
  if(dx < *xp) *xp = dx;
}
void max_double(void *max_tot,const void *test_val) {
  double dx = *(double *) test_val;
  double *xp = (double *) max_tot;
  if(dx > *xp) *xp = dx;
}

static int rd2_level = 0;
integer redistribute2_core(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			   integer nloc,integer nlocmax,bisection_data data[],
			   nodedata tree[],double c0[3],double c1[3],double volume_weight,
			   int nnei_p[1],recnlist_node **nlist_p,const double nlist_rcut,
			   const int full_tree) {
  STARTTIMER;
  static cmp_fun cmplist[] = {cmp_coord_x,cmp_coord_y,cmp_coord_z};
  const int np = (pid1-pid0)/pstride;
  const int g_pid = (pid-pid0)/pstride;
  integer pcut,nl,nr,nlt = 0 /* To keep compiler for from complaining */,nrt;
  double ccut0[3],ccut1[3];

  assert(full_tree == 0 || full_tree == 1);
  rd2_level++;

  if(np > 1) {
    double xcut;
    int idir = 0;

    {
      int i;
      for(i = 1; i<3; i++)
	if(c1[i] - c0[i] > c1[idir] - c0[idir]) idir = i;
    }

    /* Dump particle stats */ if(0) {
      double n = nloc,nmin,nmax,navg;
      MPI_Allreduce(&n,&nmin,1,MPI_DOUBLE,MPI_MIN,comm);
      MPI_Allreduce(&n,&nmax,1,MPI_DOUBLE,MPI_MAX,comm);
      MPI_Allreduce(&n,&navg,1,MPI_DOUBLE,MPI_SUM,comm);

      navg /= np << (rd2_level-1);
      if(pid == 0) {
	printf("  In %s() at %s:%d: Counts at level %d (np=%d): n = %.3e, nmin = %.3e, nmax = %.3e, navg = %.3e\n",
	       __func__,__FILE__,__LINE__,
	       rd2_level,(int) np,
	       n,nmin,nmax,navg);
	fflush(stdout);
      }
      MPI_Barrier(comm);
    }


    if(1) {
      bisection_data pivot;
      profile(RECBIS_BALANCE2,START);
      nloc = balance2(pid0,pid1,pid,pstride,comm,
		      nloc,nlocmax,sizeof(*data),data,cmplist[idir],&pivot);
      profile(RECBIS_BALANCE2,END);
      /* Determine cutting plane */ {
	double x = pivot.coord[idir],dx = fabs(x),y = x+dx,y1;
	do {
	  y1 = y;
	  dx *= 0.5;
	  y = x + dx;
	} while(y > x);
	xcut = y1;
      }
    } else {
      /*
        integer balance_hist(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
                             integer nloc,integer nlocmax,bisection_data data[],
                             integer idir,double x0,double x1,double xcut_out[1]);
      */
      const double boxVolume = (c1[0]-c0[0])*(c1[1]-c0[1])*(c1[2]-c0[2]);
      const double volume_cost = volume_weight * boxVolume;
      nloc = balance_hist(pid0,pid1,pid,pstride,comm,
			  nloc,nlocmax,data,idir,c0[idir],c1[idir],volume_cost,&xcut);
    }

    if(0) {
      int pid,np;
      MPI_Barrier(comm);
      MPI_Comm_rank(comm,&pid);
      MPI_Comm_size(comm,&np);
      if(pid == 0)
	printf("redist2_core: x0=%10.4f, xcut = %10.4f, x1=%10.4f, idir = %d\n",c0[idir],xcut,c1[idir],idir);

      {
	int nlt = 0,neq = 0,ngt = 0,nlttot,neqtot,ngttot;

	{
	  int i;
	  for(i = 0; i<nloc; i++) {
	    if(data[i].coord[idir] < xcut) nlt++;
	    if(data[i].coord[idir] == xcut) neq++;
	    if(data[i].coord[idir] > xcut) ngt++;
	  }
	}

	MPI_Allreduce(&nlt,&nlttot,1,MPI_INT,MPI_SUM,comm);
	MPI_Allreduce(&neq,&neqtot,1,MPI_INT,MPI_SUM,comm);
	MPI_Allreduce(&ngt,&ngttot,1,MPI_INT,MPI_SUM,comm);
	
	if(pid == 0) {
	  int i;
	  for(i = 0; i<np; i++) {
	    if(i > 0) {
	      MPI_Recv(&nlt,1,MPI_INT,i,39,comm,MPI_STATUS_IGNORE);
	      MPI_Recv(&neq,1,MPI_INT,i,39,comm,MPI_STATUS_IGNORE);
	      MPI_Recv(&ngt,1,MPI_INT,i,39,comm,MPI_STATUS_IGNORE);
	    }
	    printf("pid=%03d:  %8d  %8d  %8d\n",i,nlt,neq,ngt);
	  }
	  printf("*******:  %8d  %8d  %8d\n",nlttot,neqtot,ngttot);
	  
	} else {
	  MPI_Send(&nlt,1,MPI_INT,0,39,comm);
	  MPI_Send(&neq,1,MPI_INT,0,39,comm);
	  MPI_Send(&ngt,1,MPI_INT,0,39,comm);
	}
      }

      MPI_Barrier(comm);
    }

    {
      int i;
      for(i = 0; i<3; i++) {
	ccut0[i] = c0[i];
	ccut1[i] = c1[i];
      }
      ccut0[idir] = xcut;
      ccut1[idir] = xcut;
    }

    pcut = (np+1)/2;
    nl = pcut;
    nr = np - pcut;

    /* Refine neighbor list after cut... */ {
      recnlist_node sub_left,sub_right,*self,*sibling;
      int i;
      for(i = 0; i<3; i++) {
	sub_left.xx0 [i] = c0[i];
	sub_left.xx1 [i] = ccut1[i];
	sub_right.xx0[i] = ccut0[i];
	sub_right.xx1[i] = c1[i];
      }
      sub_left.pid0 = pid0;
      sub_left.pid1 = pid0+pcut*pstride;
      sub_left.stride = pstride;
      
      sub_right.pid0 = pid0+pcut*pstride;
      sub_right.pid1 = pid1;
      sub_right.stride = pstride;
      
      if(g_pid < pcut) {
	self = &sub_left;
	sibling = &sub_right;
      } else {
	self = &sub_right;
	sibling = &sub_left;
      }
      recnlist_refine(pid,comm,
		      self,sibling,
		      nlist_rcut,
		      nnei_p,nlist_p);
    }
    
    nlt = treesize(nl);
    nrt = treesize(nr);
    if(g_pid < pcut) {
      nloc = redistribute2_core(pid0,pid0+pcut*pstride,pid,pstride,comm,
				nloc,nlocmax,data,
				tree+full_tree*1,c0,ccut1,volume_weight,
				nnei_p,nlist_p,nlist_rcut,
				full_tree);
      /* Communicate subtree to/from other half of processors */ if(full_tree) {
	integer i;
	MPI_Recv(tree+1+nlt,nrt*sizeof(nodedata),MPI_BYTE,
		 pid0+(pcut+(g_pid%nr))*pstride,1023,comm,MPI_STATUS_IGNORE);      
	for(i = g_pid; i<nr; i+=nl)
	  MPI_Send(tree+1,nlt*sizeof(nodedata),MPI_BYTE,
		   pid0+(pcut+i)*pstride,1029,comm);
      }
    } else {
      nloc = redistribute2_core(pid0+pcut*pstride,pid1,pid,pstride,comm,
				nloc,nlocmax,data,
				tree+full_tree*(1+nlt),ccut0,c1,volume_weight,
				nnei_p,nlist_p,nlist_rcut,
				full_tree);
      /* Communicate subtree to/from other half of processors */ if(full_tree) {
	integer i;
	for(i = g_pid-pcut; i<nl; i+=nr)
	  MPI_Send(tree+1+nlt,nrt*sizeof(nodedata),MPI_BYTE,
		   pid0+i*pstride,1023,comm);
	MPI_Recv(tree+1,nlt*sizeof(nodedata),MPI_BYTE,
		 pid0+((g_pid-pcut)%nl)*pstride,1029,comm,MPI_STATUS_IGNORE);
      }
    }
  } else /* Single processor left... */ {

    /* Refine neighbor list until all neighbor bricks are held by one task each... */ {
      recnlist_node self;
      int i,nnei;
      for(i = 0; i<3; i++) {
	self.xx0[i] = c0[i];
	self.xx1[i] = c1[i];
      }
      self.pid0 = pid0;
      self.pid1 = pid1;
      self.stride = pstride;

      do {
	recnlist_node *nlist = nlist_p[0];
	nnei = nnei_p[0];

	for(i = 0; i<nnei; i++)
	  if((nlist[i].pid1 - nlist[i].pid0)/nlist[i].stride > 1) {
	    recnlist_refine(pid,comm,
			    &self,NULL,
			    nlist_rcut,
			    nnei_p,nlist_p);
	    break;
	  }
      } while(i < nnei);
    }

  }

  /* Create root node from sub tree data */ if(full_tree || np == 1) {
    int i;
    for(i = 0; i<3; i++) {
      tree[0].c0[i] = c0[i];
      tree[0].c1[i] = c1[i];
    }
    if(np > 1) {
      tree[0].rank = -1;
      tree[0].left = 1;
      tree[0].right = 1+nlt;
      tree[1].parent = -1;
      tree[1+nlt].parent = -(1+nlt);
    } else {
      tree[0].rank = pid;
      tree[0].left = 0;
      tree[0].right = 0;
    }
    tree[0].parent = 0;
  }

  rd2_level--;

  STOPTIMER;
  return nloc;
}


typedef struct domain_info_s {
  nodedata *tree,*neilist;
  int nnodes,self_pid,self_idx;
  int nnei;
} domain_info;

void dinf_delete(domain_info_p dip) {
  free(dip->tree);
  free(dip->neilist);
  free(dip);
}

static int overlap(double a0[3],double a1[3],double b0[3],double b1[3]) {
  int tag[3],i;
  for(i = 0; i<3; i++) {
    if(a0[i] < b0[i])
      tag[i] = (b0[i] <= a1[i]);
    else if(b0[i] < a0[i])
      tag[i] = (a0[i] <= b1[i]);
    else tag[i] = 1;
  }
  return tag[0] && tag[1] && tag[2];
}
static int inside(const double c0[3],const double c1[3],const double x[3]) {
  int tag[3],i;
  for(i = 0; i<3; i++) {
    tag[i] =
      (x[i] >= c0[i] && x[i] <= c1[i]);
  }
  return tag[0] && tag[1] && tag[2];
}
int dinf_overlap(double a0[3],double a1[3],double b0[3],double b1[3]) {
  return overlap(a0,a1,b0,b1);
}
int dinf_inside(const double c0[3],const double c1[3],const double x[3]) {
  return inside(c0,c1,x);
}

static void get_neilist_core(nodedata tree[],const nodedata self[1],double c0[3],double c1[3],
			     integer nalloc[1],integer nused[1],nodedata **list) {

  if(overlap(tree->c0,tree->c1,c0,c1)) {
    if(tree->left == 0 && tree->right == 0) {
      if(tree != self) {
	if(nused[0] >= nalloc[0]) {
	  nalloc[0] *= 2;
	  if(nused[0] >= nalloc[0]) nalloc[0] = nused[0] + 1;
	  *list = (nodedata *) realloc(*list,sizeof(nodedata) * nalloc[0]);
	  assert(*list != NULL);
	}
	(*list)[nused[0]++] = *tree;
      }
    } else {
      if(tree->left != 0)
	get_neilist_core(tree+tree->left,self,c0,c1,nalloc,nused,list);
      if(tree->right != 0)
	get_neilist_core(tree+tree->right,self,c0,c1,nalloc,nused,list);
    }
  }
}

integer dinf_get_neilist(domain_info_p dip,double rcut,nodedata **list) {
  STARTTIMER;

  assert(rcut <= getddc()->rcut &&
	 "Right now, neighbor lists are built during bisection, with a cut off equal to ddc->rcut" != NULL);

  *list = dip->neilist;

  STOPTIMER;
  return dip->nnei;
}

integer dinf_get_neilist_old(domain_info_p dip,double rcut,nodedata **list) {
  STARTTIMER;
  integer nalloc = 10,nused = 0,i;
  double c0[3],c1[3];

  assert((dip->self_idx > 0 || dip->nnodes == 1) &&
	 "Neilist only supported for full trees." != NULL);

  *list = (nodedata *) malloc(sizeof(nodedata) * nalloc);
  assert(*list != NULL);

  for(i = 0; i<3; i++) {
    c0[i] = dip->tree[dip->self_idx].c0[i] - rcut;
    c1[i] = dip->tree[dip->self_idx].c1[i] + rcut;
  }

  get_neilist_core(dip->tree,dinf_get_self(dip),c0,c1,&nalloc,&nused,list);
  *list = realloc(*list,sizeof(nodedata) * nused);
  assert(*list != NULL);

  STOPTIMER;
  return nused;
}


int nodedata_cmp_rank(const void *ap,const void *bp) {
  integer
    a = ((const nodedata *) ap)->rank,
    b = ((const nodedata *) bp)->rank;
  if(a < b) return -1;
  if(a > b) return  1;
  return 0;
}
integer dinf_get_neilist_pbc(domain_info_p dip,double rcut,nodedata **list,int pbc_flag) {
  integer nalloc = 10,nused = 0;
  double c0[3],c1[3];
  const double *box0 = dip->tree->c0,*box1 = dip->tree->c1;
  int img;

  *list = (nodedata *) malloc(sizeof(nodedata) * nalloc);
  assert(*list != NULL);

  for(img = 0; img<27; img++) {
    int d[3],dtmp = img;
    int i;
    for(i = 0; i<3; i++) {
      d[i] = (dtmp%3) - 1;
      dtmp /= 3;
      if(d[i] != 0 && ((pbc_flag>>i)&1) != 0) break;
    }
    if(i < 3) continue;

    for(i = 0; i<3; i++) {
      double imgshift = d[i]*(box1[i] - box0[i]);
      c0[i] = dip->tree[dip->self_idx].c0[i] - rcut + imgshift;
      c1[i] = dip->tree[dip->self_idx].c1[i] + rcut + imgshift;
    }
    
    get_neilist_core(dip->tree,dinf_get_self(dip),c0,c1,&nalloc,&nused,list);
    
  }

  /* Remove duplicates */ {
    nodedata *nlist = *list;
    int i,j;

    qsort(nlist,nused,sizeof(*nlist),nodedata_cmp_rank);

    j = 0;
    for(i = 1; i<nused; i++) {
      if(nlist[i].rank != nlist[j].rank) {
	j++;
	if(j < i)
	  nlist[j] = nlist[i];
      }
    }
    nused = j+1;
  }

  *list = realloc(*list,sizeof(nodedata) * nused);
  assert(*list != NULL);

  return nused;
}


const nodedata * dinf_get_self(domain_info_p dip) {
  return dip->tree + dip->self_idx;
}
static const nodedata *lookup_point_core(nodedata *tree,double x[3]) {
  assert(inside(tree->c0,tree->c1,x));

  if(tree->left == 0 && tree->right == 0)
    return tree;
  else {
    if(tree->left != 0 && inside(tree[tree->left].c0,tree[tree->left].c1,x))
      return lookup_point_core(tree + tree->left,x);
    else if(tree->right != 0 && inside(tree[tree->right].c0,tree[tree->right].c1,x))
      return lookup_point_core(tree + tree->right,x);
    else {
      assert("Children of tree not consistent, or not tesselating." == NULL);
    }
  }
}
const nodedata * dinf_lookup_point(domain_info_p dip,double x[3]) {
  assert(inside(dip->tree->c0,dip->tree->c1,x));
  return lookup_point_core(dip->tree,x);
}

static void printtree_core(nodedata *tree,int idx,int level) {
  double
    *c0 = tree[idx].c0,
    *c1 = tree[idx].c1;
  printf("%*s%3d  %4.2f,%4.2f,%4.2f\n"
	 "%*s     %4.2f,%4.2f,%4.2f\n\n",
	 20*level,"",idx,
	 c0[0],c0[1],c0[2],
	 20*level,"",
	 c1[0],c1[1],c1[2]);
  if(tree[idx].left != 0)
    printtree_core(tree,idx+tree[idx].left,level+1);
  if(tree[idx].right != 0)
    printtree_core(tree,idx+tree[idx].right,level+1);
}


void dinf_printtree(domain_info_p dip) {
  printtree_core(dip->tree,0,0);
}

#include "system.h"
extern int eam_opt_has_run;
void integer_max_accum(void *a,const void *da) {
  integer *x = (integer *) a;
  integer dx = *((const integer *) da);
  if(dx > *x) *x = dx;
}
domain_info_p redistribute2(integer np,integer pid,MPI_Comm comm,
			    integer nloc,PARTICLE pvec[],
			    double c0[3],double c1[3]) {
  if(0) if(pid == 0) {
      printf("In %s() at %s:%d: Enter\n",__func__,__FILE__,__LINE__);
      fflush(stdout);
    }

  STARTTIMER;

  static int first = 1;
  thandle t_first_call = NULL;
  if(first) {
    timestamp("Start of first bisection call.\n");
    t_first_call = rectimer_start(NULL,"redist2:first call");
  }

  {
    static thandle t_pre_barr = NULL;
    t_pre_barr = rectimer_start(t_pre_barr,"redist2:pre barr");
    MPI_Barrier(comm);
    rectimer_stop(t_pre_barr);
  }

  const int full_tree = 0;
  integer nt = full_tree ? treesize(np) : 1,nloc_r;
  nodedata *tree = (nodedata *) malloc(sizeof(nodedata) * nt);
  integer nlocmax;
  {
#if 0
    integer nmax = nloc;
    /*
    void tree_allreduce(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			integer ndata,integer sz,void *data,accum_fun acc);
    */
    tree_allreduce(0,np,pid,1,comm,
		   1,sizeof(nmax),&nmax,integer_max_accum);
#else
    int nmax,nloc_int = nloc;
    MPI_Allreduce(&nloc_int,&nmax,1,MPI_INT,MPI_MAX,comm);
#endif
    nlocmax = 3*nmax/2 + 1000;
    if(0) if(pid == 0)
      printf("In %s(): nloc=%d, nmax=%d, nlocmax=%d\n",
	     __func__,(int) nloc,(int) nmax,(int) nlocmax);
  }


  unsigned int dataBlock;
  bisection_data *data = heapGet(&dataBlock);
  //(bisection_data *) malloc(sizeof(bisection_data) * nlocmax);
  nlocmax = heapAvailable()/sizeof(*data);
  if(1) {
    int err = (data == NULL);
    int errtot;
    MPI_Allreduce(&err,&errtot,1,MPI_INT,MPI_SUM,comm);
    if(errtot > 0) {
      if(pid == 0)
	printf("At %s:%d in %s(): Out of memory, %d processes, asking for %lld slots (%lld bytes).\n",
	       __FILE__,__LINE__,__func__,errtot,(long long int) nlocmax,
	       (long long int) (sizeof(*data) * nlocmax));
      MPI_Barrier(comm);
      MPI_Finalize();
      exit(1);
    }
     
  }
  assert(data != NULL && "Not enough memory for bisection data" != NULL);

  domain_info_p dip;

  double particleWeight = 1.0,neighborWeight = 0.0,volumeWeight = 0.0;
  /* Inpect simulate object to find out weights */ 
  {
    {
      const BISECTION_BALANCE_PARMS *bisLB = (const BISECTION_BALANCE_PARMS *) (simulate_getSimulate(NULL)->ddc->loadBalance->parms);
      assert(bisLB != NULL);
      particleWeight = bisLB->particle_weight;
      neighborWeight = bisLB->neighbor_weight;
      volumeWeight   = bisLB->volume_weight;
    }
    assert(particleWeight >= 0.0 &&
	   neighborWeight >= 0.0 &&
	   volumeWeight >= 0.0 &&
	   ((particleWeight + neighborWeight + volumeWeight) > 0.0));
  }

  if(0) if(pid == 0)
    printf("In %s(): particleWeight = %.3e, neighborWeight = %.3e, g_flag = %d, eam_opt_has_run = %d\n",
	   __func__,particleWeight,neighborWeight,1,eam_opt_has_run);

  assert(tree != NULL);
  assert(data != NULL);

  /* Populate temporary coordinate array */ {
    int i;
    float *cost = system_getSystem(NULL)->collection->state->cost;
    if(0 && pid == 0) {
      if(eam_opt_has_run == 1) printf("pid=%03d @@@ Using wieght data\n",(int) pid);
      else printf("pid=%03d @@@ Setting all weights to 1\n",(int) pid);
    }
    for(i = 0; i<nloc; i++) {
      /*
        int j;
        for(j = 0; j<3; j++)
	  data[i].coord[j] = pvec[i].r[j];
      */
      data[i].coord[0] = pvec[i].r.x;
      data[i].coord[1] = pvec[i].r.y;
      data[i].coord[2] = pvec[i].r.z;
      if(0 && eam_opt_has_run == 1 && cost != NULL)
	data[i].weight = cost[i]*neighborWeight + particleWeight;
      else if(particleWeight > 0.0)
	data[i].weight = particleWeight;
      else
	data[i].weight = 1.0;
      data[i].origin = pid;
      data[i].index_at_origin = i;
      if(0) if(pid == 0 && i < 20)
	printf("i=%4d  %10.4f %10.4f %10.4f  cost = %10.4f\n",
	       i,
	       data[i].coord[0],
	       data[i].coord[1],
	       data[i].coord[2],
	       data[i].weight);
      //data[i].weight = 1.0;
    }
    eam_opt_has_run = 0;
  }

  /* Dump particle stats */ if(0) {
    double n = nloc,nmin,nmax,navg;
    MPI_Allreduce(&n,&nmin,1,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(&n,&nmax,1,MPI_DOUBLE,MPI_MAX,comm);
    MPI_Allreduce(&n,&navg,1,MPI_DOUBLE,MPI_SUM,comm);
    navg /= np;
    if(pid == 0) {
      printf("  In %s() at %s:%d: Counts on input (nloc): n = %.3e, nmin = %.3e, nmax = %.3e, navg = %.3e\n",
	     __func__,__FILE__,__LINE__,
	     n,nmin,nmax,navg);
    }
  }

  if(0) if(pid == 0) {
    printf("  In %s() at %s:%d: Entering redist2_core\n",__func__,__FILE__,__LINE__);
    fflush(stdout);
  }
  profile(RECBIS_REDIST2CORE,START);
  int nnei = 0;
  recnlist_node *nlist = NULL;
  double nlist_rcut = getddc()->rcut;

  nloc_r = redistribute2_core(0,np,pid,1,comm,
			      nloc,nlocmax,data,
			      tree,c0,c1,volumeWeight,
			      &nnei,&nlist,nlist_rcut,
			      full_tree);
  profile(RECBIS_REDIST2CORE,END);

  /* Dump neighbor list stats */ if(0) {
    double n = nnei,nmin,nmax,navg;
    MPI_Allreduce(&n,&nmin,1,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(&n,&nmax,1,MPI_DOUBLE,MPI_MAX,comm);
    MPI_Allreduce(&n,&navg,1,MPI_DOUBLE,MPI_SUM,comm);
    navg /= np;
    if(pid == 0) {
      printf("  In %s() at %s:%d: Counts after redist2 (nnei): "
	     "rcut = %.3f, n = %.3e, nmin = %.3e, nmax = %.3e, navg = %.3e\n",
	     __func__,__FILE__,__LINE__,
	     nlist_rcut,n,nmin,nmax,navg);
    }
  }


  /* Dump particle stats */ if(0) {
    double n = nloc_r,nmin,nmax,navg;
    MPI_Allreduce(&n,&nmin,1,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(&n,&nmax,1,MPI_DOUBLE,MPI_MAX,comm);
    MPI_Allreduce(&n,&navg,1,MPI_DOUBLE,MPI_SUM,comm);
    navg /= np;
    if(pid == 0) {
      printf("  In %s() at %s:%d: Counts after redist2 (nloc_r): n = %.3e, nmin = %.3e, nmax = %.3e, navg = %.3e\n",
	     __func__,__FILE__,__LINE__,
	     n,nmin,nmax,navg);
    }
  }

  if(0) if(pid == 0) {
    printf("  In %s() at %s:%d: Leave redist2_core, enter back_comm\n",__func__,__FILE__,__LINE__);
    fflush(stdout);
  }
  profile(RECBIS_BACKCOMM,START);
  back_communicate(pid,np,comm,
		   nloc_r,nlocmax,data,
		   pvec);
  profile(RECBIS_BACKCOMM,END);
  if(0) if(pid == 0) {
    printf("  In %s() at %s:%d: Leave back_comm\n",__func__,__FILE__,__LINE__);
    fflush(stdout);
  }

  //free(data);
  heapEndBlock(dataBlock,nlocmax * sizeof(*data));
  heapFree(dataBlock);

  dip = (domain_info *) malloc(sizeof(*dip));
  assert(dip != NULL);
  dip->tree = tree;
  dip->self_pid = pid;
  dip->self_idx = full_tree ? treeself(0,np,pid) : 0;
  dip->nnodes = np;

  dip->neilist = (nodedata *) malloc(sizeof(*(dip->neilist)) * nnei);
  dip->nnei = nnei;
  {
    int i;
    for(i = 0; i<nnei; i++) {
      int j;
      assert((nlist[i].pid1 - nlist[i].pid0)/nlist[i].stride == 1);
      dip->neilist[i].left = 0;
      dip->neilist[i].right = 0;
      dip->neilist[i].parent = 0;
      dip->neilist[i].rank = nlist[i].pid0;
      for(j = 0; j<3; j++) {
	dip->neilist[i].c0[j] = nlist[i].xx0[j];
	dip->neilist[i].c1[j] = nlist[i].xx1[j];
      }
    }
    free(nlist);
  }

  {
    static thandle t_barr = NULL;
    t_barr = rectimer_start(t_barr,"redist2:post barr");
    MPI_Barrier(comm);
    rectimer_stop(t_barr);
  }

  if(first) {
    rectimer_stop(t_first_call);
    first = 0;
    timestamp("End of first bisection call.\n");
  }


  STOPTIMER;

  if(0) if(pid == 0) {
    printf("In %s() at %s:%d: Leave\n",__func__,__FILE__,__LINE__);
    fflush(stdout);
  }

  return dip;
}
