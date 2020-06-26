#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#include "back_communicate.h"
#include "recnlist.h"

#include "box.h"

void recnlist_refine(integer pid,MPI_Comm comm,
		     recnlist_node *self,recnlist_node *sibling,double rcut,
		     int nnei_p[1],recnlist_node **nlist_p) {
  const int nnei = nnei_p[0];
  const double rcut_half = 0.5*rcut;
  recnlist_node
    me = *self,
    *nlist = nlist_p[0],
    *nlist2 = (recnlist_node *) malloc(sizeof(*nlist2) * (2*nnei+1));
  const integer
    pstride = self->stride,
    g_np = (self->pid1 - self->pid0) / pstride,
    g_pid = (pid - self->pid0) / pstride;
  int nnei2 = 0;

  if(sibling != NULL)
    nlist2[nnei2++] = *sibling;

  /* Communicate my data to neighbors and vice versa */ {
    int nreq = 0,nreqalloc = 6*nnei /* Two receives, at most 4 sends */;
    MPI_Request
      *reqlist = (MPI_Request *) malloc(sizeof(*reqlist) * nreqalloc);

    assert(nreqalloc == 0 || reqlist != NULL);

    /* Receive phase */ {
      int i;
      for(i = 0; i<nnei; i++) {
	/*
	  This code assumes that the neighbor represented
	  by the process group
	    inei_pid0:inei_pstride:(inei_pid1-1)
	  was split in two represented by the process groups
	    inei_pid0:inei_pstride:(inei_pid0+inei_pstride*inei_pcut-1)
	  and
	    (inei_pid0+inei_pstride*inei_pcut):inei_pstride:inei_pid1
	  this is true also for the send loop below.
	*/

	const integer
	  inei_stride = nlist[i].stride,
	  inei_pid0 = nlist[i].pid0,
	  inei_np = (nlist[i].pid1 - inei_pid0)/inei_stride,
	  inei_pcut = (inei_np + 1)/2;
	
	if(inei_np <= 1) {
	  nlist2[nnei2++] = nlist[i];
	} else {
	  MPI_Irecv(nlist2 + nnei2++,sizeof(*nlist2),MPI_BYTE,
		    inei_pid0 + inei_stride*(g_pid % inei_pcut),
		    113,comm,reqlist + nreq++);
	  MPI_Irecv(nlist2 + nnei2++,sizeof(*nlist2),MPI_BYTE,
		    inei_pid0 + inei_stride*(inei_pcut + (g_pid % (inei_np-inei_pcut))),
		    113,comm,reqlist + nreq++);
	}
      }
    }

    /* Send phase */ if(sibling != NULL) {
      int i;
      for(i = 0; i<nnei; i++) {
	const integer
	  inei_stride = nlist[i].stride,
	  inei_pid0 = nlist[i].pid0,
	  inei_np = (nlist[i].pid1 - inei_pid0)/inei_stride,
	  inei_pcut = (inei_np + 1)/2;
	integer p;
	
	for(p = g_pid; p<inei_pcut; p += g_np) {
	  assert(nreq < nreqalloc);
	  MPI_Isend(&me,sizeof(me),MPI_BYTE,
		   inei_pid0 + inei_stride*p,
		    113,comm,reqlist + nreq++);
	}
	for(p = inei_pcut+g_pid; p<inei_np; p += g_np) {
	  assert(nreq < nreqalloc);
	  MPI_Isend(&me,sizeof(me),MPI_BYTE,
		   inei_pid0 + inei_stride*p,
		   113,comm,reqlist + nreq++);
	}
      }
    }

    MPI_Waitall(nreq,reqlist,MPI_STATUSES_IGNORE);
    free(reqlist);
  }

  /* Prune new neighbor list */ {
    THREE_VECTOR c0v = box_get_corner(NULL),c1v = box_get_urcorner(NULL);
    double boxl[3] = { c1v.x-c0v.x , c1v.y-c0v.y , c1v.z-c0v.z };
    int pbc = box_get_boundary_conditions(NULL);

    assert(box_getType(NULL) == ORTHORHOMBIC);		 
    assert(pbc >= 0 && pbc <= 7);

    recnlist_node mebig = me;
    int i,j = 0;
    
    /* A self copy with skin added */ {
      int k;
      for(k = 0; k<3; k++) {
	mebig.xx0[k] -= rcut_half;
	mebig.xx1[k] += rcut_half;
      }
    }
    
    for(i = 0; i<nnei2; i++) {
      int is_nei = 0;
      recnlist_node
	tmpbig = nlist2[i];

      /* Add skin to neighbor domain */ {
	int k;
	for(k = 0; k<3; k++) {
	  tmpbig.xx0[k] -= rcut_half;
	  tmpbig.xx1[k] += rcut_half;
	}
      }

      /* Could also symmetrize using pid0 nlist[2].pid0... */

      if(dinf_overlap(mebig.xx0,mebig.xx1,tmpbig.xx0,tmpbig.xx1))
	is_nei = 1;
      else if(pbc) {
	recnlist_node tmpimg = tmpbig,meimg = mebig;
	const int
	  xlim = (pbc>>0)&1,
	  ylim = (pbc>>1)&1,
	  zlim = (pbc>>2)&1;
	int ix,iy,iz;
	for(ix = -xlim; ix<=xlim; ix++) {
	  const double dx = 0.5*boxl[0]*ix;
	  tmpimg.xx0[0] = tmpbig.xx0[0] + dx;
	  tmpimg.xx1[0] = tmpbig.xx1[0] + dx;
	  meimg.xx0 [0] = mebig.xx0 [0] - dx;
	  meimg.xx1 [0] = mebig.xx1 [0] - dx;
	  for(iy = -ylim; iy<=ylim; iy++) {
	    const double dy = 0.5*boxl[1]*iy;
	    tmpimg.xx0[1] = tmpbig.xx0[1] + dy;
	    tmpimg.xx1[1] = tmpbig.xx1[1] + dy;
	    meimg.xx0 [1] = mebig.xx0 [1] - dy;
	    meimg.xx1 [1] = mebig.xx1 [1] - dy;

	    for(iz = -zlim; iz<=zlim; iz++) {
	      const double dz = 0.5*boxl[2]*iz;
	      tmpimg.xx0[2] = tmpbig.xx0[2] + dz;
	      tmpimg.xx1[2] = tmpbig.xx1[2] + dz;
	      meimg.xx0 [2] = mebig.xx0 [2] - dz;
	      meimg.xx1 [2] = mebig.xx1 [2] - dz;
	      
	      if(dinf_overlap(meimg.xx0,meimg.xx1,tmpimg.xx0,tmpimg.xx1)) {
		is_nei = 1;
		goto out;
	      }

	    }
	  }
	}
      }
    out:
      if(is_nei) {
	if(i > j)
	  nlist2[j] = nlist2[i];
	j++;
      }
    }
    nnei2 = j;
  }

  free(nlist);
  nlist2 = (recnlist_node *) realloc(nlist2,sizeof(*nlist2) * nnei2);
  nlist_p[0] = nlist2;
  nnei_p[0] = nnei2;
}
