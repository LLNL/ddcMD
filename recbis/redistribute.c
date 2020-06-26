/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <stdio.h>

#include "treecomm.h"

#include "redistribute.h"
#include "balance.h"
/*
integer balance(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		integer nloc,integer nlocmax,integer sz,void *data,cmp_fun cmp);
*/


static int level = 0;
integer redistribute(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc,integer nlocmax,integer sz,void *data,
		     int ndir,int idir,cmp_fun cmp_coord[/*ndir*/]) {
  const int debug = 0,serialize = debug;
  const int np = (pid1-pid0)/pstride;
  const int g_pid = (pid-pid0)/pstride;
  integer pcut;

  if(np <= 1) return nloc;

  level+=1;

  nloc = balance(pid0,pid1,pid,pstride,comm,
		 nloc,nlocmax,sz,data,cmp_coord[idir]);

  if(debug) {
    integer ntot = nloc,i;
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(ntot),&ntot,sum_integer);
    if(g_pid == 0)
      printf("1%03d%04d1 %*s Balance in redistribute(): ntot = %3d  np = %2d  pid = %2d\n\n",
	     (int) np,(int) pid0,2*level,"",
	     (int) ntot,(int) np,(int) pid);
    if(0) for(i = 0; i<np; i++) {
      if(i == g_pid) {
	printf("1%03d%04d2 %*s   pid=%2d  g_pid=%2d  nloc=%d\n",
	       (int) np,(int) pid0,2*level,"",
	       (int) pid,(int) g_pid,(int) nloc);
      }
      /* Barrier... */
      tree_allreduce(pid0,pid1,pid,pstride,comm,
		     1,sizeof(ntot),&ntot,sum_integer);
    }
  }
  
  pcut = (np+1)/2;
  if(g_pid < pcut) {
    nloc = redistribute(pid0,pid0+pcut*pstride,pid,pstride,comm,
			nloc,nlocmax,sz,data,ndir,(idir+1)%ndir,cmp_coord);
    if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);
  } else {
    if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);
    nloc = redistribute(pid0+pcut*pstride,pid1,pid,pstride,comm,
			nloc,nlocmax,sz,data,ndir,(idir+1)%ndir,cmp_coord);
  }

  level-=1;
  if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);

  return nloc;
}

