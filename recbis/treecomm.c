/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/


#include <stdlib.h>

#include <mpi.h>
#include "treecomm.h"

#include "genptr.h"

#include "rectimer.h"
#include "functimer.h"


integer tree_gather_inc(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer nloc_max,integer nloc,integer sz,void *data) {
  integer gsize = (pid1-pid0)/pstride;
  integer g_pid = (pid-pid0)/pstride;
  integer bit,nacc;

  nacc = nloc;
  for(bit = 1;  bit<gsize; bit*=2)
    if((g_pid & bit) > 0) {
      MPI_Send(data,nacc*sz,MPI_BYTE,pid-bit*pstride,993,comm);
      break;
    } else if(g_pid+bit < gsize) {
      MPI_Status status;
      integer nget = gsize - (g_pid + bit);
      if(nget > bit) nget = bit;
      nget *= nloc_max;
      MPI_Recv(idx_ptr(data,nacc,sz),nget*sz,MPI_BYTE,pid+bit*pstride,993,comm,&status);
      /* nget not be an int, so use a temporary for the call to MPI_Get_count() */ {
	int m;
	MPI_Get_count(&status,MPI_BYTE,&m);
	nget = m/sz;
      }
      nacc += nget;
    }
  return nacc;
}


void tree_gather(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc,integer sz,void *data) {
  STARTTIMER;
  integer gsize = (pid1-pid0)/pstride;
  integer g_pid = (pid-pid0)/pstride;
  integer bit,nacc;

  nacc = nloc;
  for(bit = 1;  bit<gsize; bit*=2)
    if((g_pid & bit) > 0) {
      MPI_Send(data,nacc*sz,MPI_BYTE,pid-bit*pstride,993,comm);
      break;
    } else if(g_pid+bit < gsize) {
      MPI_Status status;
      integer nget = gsize - (g_pid + bit);
      if(nget > bit) nget = bit;
      nget *= nloc;
      MPI_Recv(idx_ptr(data,nacc,sz),nget*sz,MPI_BYTE,pid+bit*pstride,993,comm,&status);
      nacc += nget;
    }
  STOPTIMER;
}

void tree_broadcast(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer ndata,integer sz,void *data) {
  STARTTIMER;
  integer gsize = (pid1-pid0)/pstride;
  integer g_pid = (pid-pid0)/pstride;

  if(g_pid > 0) {
    MPI_Status status;
    MPI_Recv(data,ndata*sz,MPI_BYTE,pid0 + ((g_pid-1)/2)*pstride,211,comm,&status);
  }
  if(2*g_pid+1 < gsize)
    MPI_Send(data,ndata*sz,MPI_BYTE,pid0+(2*g_pid+1)*pstride,211,comm);
  if(2*g_pid+2 < gsize)
    MPI_Send(data,ndata*sz,MPI_BYTE,pid0+(2*g_pid+2)*pstride,211,comm);
  STOPTIMER;
}

void tree_reduce(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer ndata,integer sz,void *data,accum_fun acc) {
  STARTTIMER;
  integer gsize = (pid1-pid0)/pstride;
  integer g_pid = (pid-pid0)/pstride;
  integer bit;

  void *temp = (void *) malloc(sz * ndata);

  for(bit = 1;  bit<gsize; bit*=2) {
    MPI_Status status;
    integer j;

    if((g_pid & bit) > 0) {
      MPI_Send(data,ndata*sz,MPI_BYTE,pid-bit*pstride,994,comm);
      break;
    } else if(g_pid + bit < gsize) {
      MPI_Recv(temp,ndata*sz,MPI_BYTE,pid+bit*pstride,994,comm,&status);
      for(j = 0; j<ndata; j++) {
	//data[j] = comm_func(data[j],temp[j]);
	acc(idx_ptr(data,j,sz),idx_ptr(temp,j,sz));
      }
    }

  }

  free(temp);
  STOPTIMER;
}

void tree_allreduce(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer ndata,integer sz,void *data,accum_fun acc) {
  STARTTIMER;
  tree_reduce(pid0,pid1,pid,pstride,comm,
	      ndata,sz,data,acc);
  tree_broadcast(pid0,pid1,pid,pstride,comm,
		 ndata,sz,data);
  STOPTIMER;
}

void tree_barrier(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm) {
  STARTTIMER;
  integer x = 1;
  tree_allreduce(pid0,pid1,pid,pstride,comm,
		 1,sizeof(x),&x,sum_integer);
  if(x*pstride != pid1-pid0) {
    MPI_Abort(comm,7);
  }
  STOPTIMER;
}
