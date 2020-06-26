/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "exchange.h"

#include "genptr.h"

#include "txmap.h"

void partition(integer n,integer np,integer pid,
	       integer idx0[1],integer idx1[1]) {
  integer q = n/np,r = n%np;
  if(np <= 0) {
    printf("Error in call to partition(): n=%d, np=%d, pid=%d\n",
	   (int) n,(int) np,(int) pid);
  }
  assert(np > 0);

  if(pid < r) {
    idx0[0] = q*pid + pid;
    idx1[0] = idx0[0] + q + 1;
  } else {
    idx0[0] = q*pid + r;
    idx1[0] = idx0[0] + q;
  }
}

void exchange_sched(integer pid0,integer pid1,integer pid,integer pstride,
		    integer rp0,integer rp1,integer rpstride,
		    MPI_Comm comm,
		    integer nhavetot,integer ngoaltot,
		    integer rnhavetot,integer rngoaltot,
		    txmap_struct *txmap) {
  const integer np = (pid1-pid0)/pstride,rnp = (rp1-rp0)/rpstride;
  const integer g_pid = (pid-pid0)/pstride;

  integer nloc;
  integer idx_have_0,idx_have_1,idx_goal_0,idx_goal_1;
  integer idx_trans_0,idx_trans_1;
  integer clist[4],nlist[3];

  if(rnp == 0) return;

  partition(nhavetot,np,g_pid,&idx_have_0,&idx_have_1);
  partition(ngoaltot,np,g_pid,&idx_goal_0,&idx_goal_1);

  nloc = idx_have_1 - idx_have_0;

  idx_trans_0 = abs_integer(idx_have_0 - idx_goal_0);
  idx_trans_1 = abs_integer(idx_have_1 - idx_goal_1);

  /* Calculate break points for changes in the number of transmitted items. */ {
    const integer q1 = rngoaltot/rnp,r1 = rngoaltot%rnp;
    const integer q2 = rnhavetot/rnp,r2 = rnhavetot%rnp;
    integer n0a,nab,nb1,ca,cb;

    n0a = abs_integer( (q1+1) - (q2+1) );
    nb1 = abs_integer( q1 - q2 );
    
    if(r1 < r2) {
      ca = r1;
      cb = r2;
      nab = abs_integer( q1 - (q2+1) );
    } else {
      ca = r2;
      cb = r1;
      nab = abs_integer( (q1+1) - q2 );
    }
    
    clist[0] = 0;  nlist[0] = n0a;
    clist[1] = ca; nlist[1] = nab;
    clist[2] = cb; nlist[2] = nb1;
    clist[3] = rnp;
  }

  /* Communicate items */ {
    integer noff = 0,kk;
    for(kk = 0; kk<3; kk++) {
      if(nlist[kk] == 0) continue;
      while(idx_trans_0 < idx_trans_1 && idx_trans_0 < noff + (clist[kk+1]-clist[kk])*nlist[kk]) {
	const integer g_rp = clist[kk] + (idx_trans_0-noff)/nlist[kk];
	const integer ntrans =
	  min_integer(idx_trans_1 - idx_trans_0,
		      (nlist[kk]*(g_rp + 1 - clist[kk]) + noff) -
		      idx_trans_0);
	
	if(ntrans > 0) {
	  if(txmap->nused >= txmap->nalloc) {
	    txmap->nalloc = (txmap->nalloc * 3)/2 + 10;
	    txmap->list = (txmap_entry *)
	      realloc(txmap->list,sizeof(*(txmap->list)) * txmap->nalloc);
	  }
	  
	  if(nhavetot > ngoaltot) {
	    txmap->list[txmap->nused  ].count = ntrans;
	    txmap->list[txmap->nused++].pid   = rp0 + g_rp*rpstride;
	    /*
	    printf("%d sending [ %d , %d ] = %d items to %d\n",
		   (int) pid,
		   (int) (nloc-ntrans),
		   (int) (nloc-1),
		   (int) ntrans,
		   (int) (rp0+g_rp*rpstride));
	    */
	    nloc -= ntrans;
	  } else {
	    txmap->list[txmap->nused  ].count = -ntrans;
	    txmap->list[txmap->nused++].pid   = rp0 + g_rp*rpstride;
	    /*
	    printf("%d receiving [ %d , %d ] = %d items from %d\n",
		   (int) pid,
		   (int) nloc,
		   (int) (nloc+ntrans-1),
		   (int) ntrans,
		   (int) (rp0+g_rp*rpstride));
	    */
	    nloc += ntrans;
	  }
	}
	idx_trans_0 += ntrans;
	
      }
      noff += (clist[kk+1] - clist[kk])*nlist[kk];
    }
  }

  if(nloc != idx_goal_1 - idx_goal_0) {
    printf("{%s:%d} pid=%03d: nloc=%d idx_goal=[%d-%d]=%d\n",
	   __FILE__,__LINE__,
	   (int) pid,
	   (int) nloc,
	   (int) idx_goal_0,
	   (int) (idx_goal_1-1),
	   (int) (idx_goal_1-idx_goal_0));
    MPI_Abort(comm,19);
  }
  //nloc_p[0] = nloc;
}

void exchange(integer pid0,integer pid1,integer pid,integer pstride,
	      integer rp0,integer rp1,integer rpstride,
	      MPI_Comm comm,
	      integer nhavetot,integer ngoaltot,
	      integer rnhavetot,integer rngoaltot,
	      integer nlocmax,integer sz,void *data) {
  const integer np = (pid1-pid0)/pstride,rnp = (rp1-rp0)/rpstride;
  const integer g_pid = (pid-pid0)/pstride;

  integer nloc;
  integer idx_have_0,idx_have_1,idx_goal_0,idx_goal_1;
  integer idx_trans_0,idx_trans_1;
  integer clist[4],nlist[3];

  if(rnp == 0) return;

  partition(nhavetot,np,g_pid,&idx_have_0,&idx_have_1);
  partition(ngoaltot,np,g_pid,&idx_goal_0,&idx_goal_1);

  if(idx_goal_1 - idx_goal_0 > nlocmax) {
    printf("{%s:%d} pid=%03d: nlocgoal=%d nlocmax=%d\n",
	   __FILE__,__LINE__,
	   (int) pid,
	   (int) (idx_goal_1 - idx_goal_0),
	   (int) nlocmax);
    MPI_Abort(comm,19);
  }

  nloc = idx_have_1 - idx_have_0;

  idx_trans_0 = abs_integer(idx_have_0 - idx_goal_0);
  idx_trans_1 = abs_integer(idx_have_1 - idx_goal_1);

  /* Calculate break points for changes in the number of transmitted items. */ {
    const integer q1 = rngoaltot/rnp,r1 = rngoaltot%rnp;
    const integer q2 = rnhavetot/rnp,r2 = rnhavetot%rnp;
    integer n0a,nab,nb1,ca,cb;

    n0a = abs_integer( (q1+1) - (q2+1) );
    nb1 = abs_integer( q1 - q2 );
    
    if(r1 < r2) {
      ca = r1;
      cb = r2;
      nab = abs_integer( q1 - (q2+1) );
    } else {
      ca = r2;
      cb = r1;
      nab = abs_integer( (q1+1) - q2 );
    }
    
    clist[0] = 0;  nlist[0] = n0a;
    clist[1] = ca; nlist[1] = nab;
    clist[2] = cb; nlist[2] = nb1;
    clist[3] = rnp;
  }

  /* Communicate items */ {
    integer noff = 0,kk;
    for(kk = 0; kk<3; kk++) {
      if(nlist[kk] == 0) continue;
      while(idx_trans_0 < idx_trans_1 && idx_trans_0 < noff + (clist[kk+1]-clist[kk])*nlist[kk]) {
	const integer g_rp = clist[kk] + (idx_trans_0-noff)/nlist[kk];
	const integer ntrans =
	  min_integer(idx_trans_1 - idx_trans_0,
		      (nlist[kk]*(g_rp + 1 - clist[kk]) + noff) -
		      idx_trans_0);
	
	if(ntrans > 0) {
	  if(nhavetot > ngoaltot) {
	    MPI_Send(idx_ptr(data,nloc-ntrans,sz),ntrans*sz,MPI_BYTE,rp0+g_rp*rpstride,719,comm);
	    /*
	    printf("%d sending [ %d , %d ] = %d items to %d\n",
		   (int) pid,
		   (int) (nloc-ntrans),
		   (int) (nloc-1),
		   (int) ntrans,
		   (int) (rp0+g_rp*rpstride));
	    */
	    nloc -= ntrans;
	  } else {
	    MPI_Status status;
	    MPI_Recv(idx_ptr(data,nloc,sz),ntrans*sz,MPI_BYTE,rp0+g_rp*rpstride,719,comm,&status);
	    /*
	    printf("%d receiving [ %d , %d ] = %d items from %d\n",
		   (int) pid,
		   (int) nloc,
		   (int) (nloc+ntrans-1),
		   (int) ntrans,
		   (int) (rp0+g_rp*rpstride));
	    */
	    nloc += ntrans;
	  }
	}
	idx_trans_0 += ntrans;
	
      }
      noff += (clist[kk+1] - clist[kk])*nlist[kk];
    }
  }

  if(nloc != idx_goal_1 - idx_goal_0) {
    printf("{%s:%d} pid=%03d: nloc=%d idx_goal=[%d-%d]=%d\n",
	   __FILE__,__LINE__,
	   (int) pid,
	   (int) nloc,
	   (int) idx_goal_0,
	   (int) (idx_goal_1-1),
	   (int) (idx_goal_1-idx_goal_0));
    MPI_Abort(comm,19);
  }
  //nloc_p[0] = nloc;
}
