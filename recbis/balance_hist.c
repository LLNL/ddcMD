/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>

#include "genptr.h"
#include "equalize.h"
#include "exchange.h"

#include "bisection_data.h"
#include "approx_median.h"

#include "balance_hist.h"

#include "rectimer.h"

static void sum_max_min_integer(void *va,const void *vda) {
  integer *a = (integer *) va;
  const integer *da = (const integer *) vda;

  a[0] += da[0];
  if(da[1] > a[1]) a[1] = da[1];
  if(da[2] < a[2]) a[2] = da[2];
}

integer balance_hist(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc,integer nlocmax,bisection_data data[],
		     integer idir,double x0,double x1,double volume_cost,
		     double xcut_out[1]) {
  static thandle t_all = NULL;
  const integer sz = sizeof(*data);

  const integer np = (pid1 - pid0)/pstride;
  const integer g_pid = (pid - pid0)/pstride;
  const integer pcut = (np+1)/2;
  integer ntot,nkeep,nsend,nrecv;

  const int glob_barr = 0;

  t_all = rectimer_start(t_all,__func__);

  if(np == 1) {
    xcut_out[0] = x1;
    rectimer_stop(t_all);
    return nloc;
  }

  //assert((np & (np-1)) == 0);
  if(glob_barr) {
    static thandle barr0 = NULL;
    barr0 = rectimer_start(barr0,"balhist:barr0");
    MPI_Barrier(comm);
    rectimer_stop(barr0);
  }

  //printf("pid=%03d: Entry.\n",(int) pid);

  /* Count total number of particles, and figure out whether we need to run equalize() */ {
    integer sum_max_min[3],maxnloc,minbufspace;
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
      xcut_out[0] = x1;
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

  if(glob_barr) {
    static thandle barr1 = NULL;
    barr1 = rectimer_start(barr1,"balhist:barr1");
    MPI_Barrier(comm);
    rectimer_stop(barr1);
  }

  //printf("pid=%03d: Call median finder.\n",(int) pid);

  /* Find weighted median, and partition data */ {
    const double frac = ((double) pcut)/(double) np;
    const double xcut = approx_kstat(pid0,pid1,pid,pstride,comm,
				     nloc,data,idir,x0,x1,volume_cost,frac);
    xcut_out[0] = xcut;
    
    /* Local partitioning according to median (xcut) */ {
      integer i = 0,j = nloc-1;

      if(g_pid < pcut) {
	while(i <= j) {
	  if(data[i].coord[idir] <= xcut) i++;
	  else if(data[j].coord[idir] > xcut) j--;
	  else {
	    bisection_data tmp = data[i];
	    data[i] = data[j];
	    data[j] = tmp;
	    i++;
	    j--;
	  }
	}
      } else {
	while(i <= j) {
	  if(! (data[i].coord[idir] <= xcut)) i++;
	  else if(! (data[j].coord[idir] > xcut)) j--;
	  else {
	    bisection_data tmp = data[i];
	    data[i] = data[j];
	    data[j] = tmp;
	    i++;
	    j--;
	  }
	}
      }
      nkeep = i;
      nsend = nloc - nkeep;
    }
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

  if(glob_barr) {
    static thandle barr4 = NULL;
    barr4 = rectimer_start(barr4,"balhist:barr4");
    MPI_Barrier(comm);
    rectimer_stop(barr4);
  }

  //printf("pid=%03d: Leave.\n",(int) pid);

  rectimer_stop(t_all);
  return nkeep + nrecv;
}

#ifdef COMPILE_BALANCE_TEST_MAIN
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[]) {
  int pid,np,nmax,n;
  bisection_data *v;
  double cut;
  int i;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);

  srand(pid*3331);
  nmax = atoi(argv[1]);
  n = rand() % nmax;

  {
    int ntot;
    MPI_Allreduce(&n,&ntot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(pid == 0)
      printf("nmax = %d (%d), ntot = %d (%.3f %%)\n",nmax,nmax*np,ntot,((double) ntot)/(nmax*np));
    MPI_Barrier(MPI_COMM_WORLD);
  }

  v = (bisection_data *) malloc(sizeof(*v) * 2*nmax);

  for(i = 0; i<n; i++) {
    v[i].coord[0] = rand() / (double) RAND_MAX;
    v[i].weight = 1.0;
    v[i].origin = pid;
    v[i].index_at_origin = i;
  }
  /*
  integer balance_hist(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		       integer nloc,integer nlocmax,bisection_data data[],
		       integer idir,double x0,double x1,double xcut_out[1]);
  */
  n = balance_hist(0,np,pid,1,MPI_COMM_WORLD,
		   n,2*nmax,v,0,0.0,1.0,&cut);

  {
    int nlt = 0,neq = 0,ngt = 0,nlttot,neqtot,ngttot;
    for(i = 0; i<n; i++) {
      if(v[i].coord[0] < cut) nlt++;
      if(v[i].coord[0] == cut) neq++;
      if(v[i].coord[0] > cut) ngt++;
    }
    
    MPI_Allreduce(&nlt,&nlttot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&neq,&neqtot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&ngt,&ngttot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(pid == 0) {
      for(i = 0; i<np; i++) {
	if(i > 0) {
	  MPI_Recv(&nlt,1,MPI_INT,i,39,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Recv(&neq,1,MPI_INT,i,39,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Recv(&ngt,1,MPI_INT,i,39,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	printf("pid=%03d:  %8d  %8d  %8d\n",i,nlt,neq,ngt);
      }
      printf("*******:  %8d  %8d  %8d\n",nlttot,neqtot,ngttot);

    } else {
      MPI_Send(&nlt,1,MPI_INT,0,39,MPI_COMM_WORLD);
      MPI_Send(&neq,1,MPI_INT,0,39,MPI_COMM_WORLD);
      MPI_Send(&ngt,1,MPI_INT,0,39,MPI_COMM_WORLD);
    }
  }

  free(v);

  MPI_Barrier(MPI_COMM_WORLD);
  if(pid == 0)
    printf("Finished.\n");
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;

}

#endif
