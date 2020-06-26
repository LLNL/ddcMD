/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <mpi.h>


#include "genptr.h"
#include "integer.h"
#include "treecomm.h"

#include "equalize.h"
#include "exchange.h"
#include "balance.h"

#include "rectimer.h"
#include "functimer.h"

#include "alltoall_sparse.h"

static void sum_max_min_integer(void *va,const void *vda) {
  integer *a = (integer *) va;
  const integer *da = (const integer *) vda;

  a[0] += da[0];
  if(da[1] > a[1]) a[1] = da[1];
  if(da[2] < a[2]) a[2] = da[2];
}

integer alltoall_sparse(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			integer nloc,integer nlocmax,
			integer sz,void *data,integer (*get_dest)(const void *element)) {
  STARTTIMER;
  const integer np = (pid1 - pid0)/pstride;
  const integer g_pid = (pid - pid0)/pstride;
  const integer pcut = (np + 1)/2;
  integer nkeep,nsend,nrecv;

  if(np <= 1) { STOPTIMER; return nloc; }

  //printf("pid=%03d: Entry.\n",(int) pid);

  /* Count total number of particles, and figure out whether we need to run equalize() */ {
    integer sum_max_min[3],ntot,maxnloc,minbufspace;
    sum_max_min[0] = nloc;
    sum_max_min[1] = nloc;
    sum_max_min[2] = nlocmax - nloc;

    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(sum_max_min),sum_max_min,sum_max_min_integer);
    ntot = sum_max_min[0];
    maxnloc = sum_max_min[1];
    minbufspace = sum_max_min[2];

    if(ntot == 0) {
      /* Early bilout in case there are no particles */
      rectimer_stop(t_all);
      return nloc;
    }

    //printf("pid=%03d: ntot=%lld, mxnloc=%lld, minbufspace=%lld.\n",(int) pid,
    //	   ntot,maxnloc,minbufspace);

    if(maxnloc*np > 2*ntot || minbufspace <= 2*maxnloc) {
      /* Run equalize() if the data amount imbalance is
	 too big or some node is short on buffer space */

      if(minbufspace < maxnloc) {
	/* This is not a definitive test, I don't have a proof yet for a slim upper bound
	   on sufficient buffer size to gurantee that dead-lock is impossible in equalize().
	   I could call the old (slow) code in this case, since it has reqsonable known
	   bounds. */
	if(g_pid == 0) {
	  time_t now = time(NULL);
	  printf("pid=%8lld: Time is %s. In %s() at %s:%d: Warning: Low on buffer space, "
		 "dead-lock possible. If you don't see a message that equalize() succeded "
		 "within a minute, please rerun with bigger heap.\n",
		 (long long int) pid,ctime(&now),__func__,__FILE__,__LINE__);
	  fflush(stdout);
	}
      }

      //printf("pid=%03d: Calling equalize.\n",(int) pid);
      /* Call equalize() */ {
	const integer ntot2 = equalize(pid0,pid1,pid,pstride,comm,
				       &nloc,nlocmax,sz,data);
	assert(ntot2 == ntot);
      }
      //printf("pid=%03d: Equalize returned.\n",(int) pid);

      if(minbufspace < maxnloc) {
	tree_barrier(pid0,pid1,pid,pstride,comm);
	if(g_pid == 0) {
	  time_t now = time(NULL);
	  printf("pid=%8lld: Time is %s. In %s() at %s:%d: Warning clear: equalize() succeeded!\n",
		 (long long int) pid,ctime(&now),__func__,__FILE__,__LINE__);
	  fflush(stdout);
	}
      }
    }
  }


  /* Put elements less than pcut first, the rest last */ {
    const integer dest_cut = pid0 + pcut*pstride;
    integer i = 0,j = nloc - 1;

    if(g_pid < pcut) {
      while(i <= j) {
	if(get_dest(idx_ptr(data,i,sz)) < dest_cut) i++;
	else if(get_dest(idx_ptr(data,j,sz)) >= dest_cut) j--;
	else swap(data,i,j,sz);
      }
    } else {
      while(i <= j) {
	if(! (get_dest(idx_ptr(data,i,sz)) < dest_cut)) i++;
	else if(! (get_dest(idx_ptr(data,j,sz)) >= dest_cut)) j--;
	else swap(data,i,j,sz);
      }
    }
    nkeep = i;
    nsend = nloc - nkeep;
  }


  //printf("pid=%03d: Data exchange.\n",(int) pid);

  /* Exchange data */ {
    if(g_pid < pcut) {
      MPI_Request reqlist[3];
      int nreq = 0;

      integer rpid = pcut + g_pid;
      if(rpid >= np) rpid--;
      rpid = pid0 + pstride*rpid;
      MPI_Isend(&nsend,sizeof(nsend),MPI_BYTE,rpid,
		343,comm,reqlist+nreq++);
      if(g_pid < np-pcut) {
	MPI_Recv(&nrecv,sizeof(nrecv),MPI_BYTE,rpid,
		  343,comm,MPI_STATUS_IGNORE);
	//printf("pid=%03d: recv %3d from %3d.\n",(int) pid,(int) nrecv,(int) rpid);
        assert(nrecv <= nlocmax - nloc);
	if(nrecv > 0)
	  MPI_Irecv(idx_ptr(data,nloc,sz),nrecv*sz,MPI_BYTE,rpid,
		    342,comm,reqlist+nreq++);
      } else
	nrecv = 0;
      //printf("pid=%03d: send %3d  to  %3d.\n",(int) pid,(int) nsend,(int) rpid);
      if(nsend > 0)
	MPI_Isend(idx_ptr(data,nkeep,sz),nsend*sz,MPI_BYTE,rpid,
		  342,comm,reqlist+nreq++);
      //printf("pid=%03d: Waitall...\n",(int) pid);
      MPI_Waitall(nreq,reqlist,MPI_STATUSES_IGNORE);
    } else {
      MPI_Request reqlist[4],rreqlist[2];
      int nreq = 0,nrreq = 0;
      integer rcvec[2];

      const integer
	g_rpid = g_pid - pcut,
	rpid = pid0 + pstride*g_rpid;

      MPI_Isend(&nsend,sizeof(nsend),MPI_BYTE,rpid,
		343,comm,reqlist+nreq++);
      MPI_Irecv(rcvec+0,sizeof(*rcvec),MPI_BYTE,rpid,
		343,comm,rreqlist+nrreq++);
      if(g_pid == np-1 && g_rpid+1<pcut)
	MPI_Irecv(rcvec+1,sizeof(*rcvec),MPI_BYTE,
		  rpid+1,343,comm,rreqlist+nrreq++);

      /* Receive data */ {
	nrecv = 0;
	int i;
	for(i = 0; i<nrreq; i++) {
	  MPI_Wait(rreqlist+i,MPI_STATUS_IGNORE);
          assert(nrecv+rcvec[i] <= nlocmax - nloc);
	  //printf("pid=%03d: recv %3d from %3d.\n",(int) pid,(int) rcvec[i],(int) g_rpid+i);
	  if(rcvec[i] > 0)
	    MPI_Irecv(idx_ptr(data,nloc+nrecv,sz),rcvec[i]*sz,MPI_BYTE,pid0 + pstride*(g_rpid+i),
		      342,comm,reqlist+nreq++);
	  nrecv += rcvec[i];
	}
      }
      //printf("pid=%03d: send %3d  to  %3d.\n",(int) pid,(int) nsend,(int) rpid);
      if(nsend > 0)
	MPI_Isend(idx_ptr(data,nkeep,sz),nsend*sz,MPI_BYTE,rpid,
		  342,comm,reqlist+nreq++);
      //printf("pid=%03d: Waitall...\n",(int) pid);
      MPI_Waitall(nreq,reqlist,MPI_STATUSES_IGNORE);
    }

    //printf("pid=%03d: Compact local data\n",(int) pid);
    /* Compact local data by copying received data into hole created by sends */ {
      if(nrecv < nsend)
	memcpy(idx_ptr(data,nkeep,sz),
	       idx_ptr(data,nloc,sz),
	       nrecv * sz);
      else
	memcpy(idx_ptr(data,nkeep,sz),
	       idx_ptr(data,nloc+nrecv-nsend,sz),
	       nsend * sz);
    }
  }

  nloc = nkeep + nrecv;

  STOPTIMER;
  /* Return number of local elements (after recursion) */ {
    if(g_pid < pcut) {
      return alltoall_sparse(pid0,pid0+pcut*pstride,pid,pstride,comm,
			       nloc,nlocmax,sz,data,get_dest);
    } else {
      return alltoall_sparse(pid0+pcut*pstride,pid1,pid,pstride,comm,
			     nloc,nlocmax,sz,data,get_dest);
    }
  }
}
