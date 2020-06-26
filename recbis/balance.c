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

#include "balance.h"

#include "treecomm.h"
#include "parkstat.h"
#include "equalize.h"
#include "exchange.h"

#include "rectimer.h"
#include "functimer.h"


static void vprint(int pid,int n1,int n2,double *data) {
  int i;

  printf("pid=%03d: (%d)",pid,n1);
  for(i = 0; i<n1; i++)
    printf("  %.2f",data[i]);
  printf(" ## (%d) ",n2);
  for(i = 0; i<n2; i++)
    printf("  %.2f",data[n1+i]);
  printf("\n");
}

integer balance2(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc,integer nlocmax,integer sz,void *data,cmp_fun cmp, void *pivot) {

  const int debug = 0,serialize = 0;
  const integer np = (pid1 - pid0)/pstride;
  const integer g_pid = (pid - pid0)/pstride;
  integer ntot,nle,ngt,ngt_tot,nle_tot,npivots;
  integer pcut,kcut;
  /*void *pivot;*/

  if(np == 1) return nloc;

  STARTTIMER;

  ntot = nloc;
  tree_allreduce(pid0,pid1,pid,pstride,comm,
		 1,sizeof(integer),&ntot,sum_integer);

  if(ntot == 0) { STOPTIMER; return nloc; }
  pcut = (np+1)/2;
  kcut = (ntot*pcut-(np+1)/2)/np;


  if(debug) {
    integer x = 1;
    /* Barrier... */
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(x),&x,sum_integer);
    if(g_pid == 0)
      printf("{%s:%d} pid=%03d: nloc=%d ntot=%d np=%d pcut=%d kcut=%d, calling parkstat\n",
	     __FILE__,__LINE__,
	     (int) pid,(int) nloc,(int) ntot,(int) np,(int) pcut,(int) kcut);
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(x),&x,sum_integer);
  }

  /*pivot = malloc(sz);*/
  npivots = parkstat(pid0,pid1,pid,pstride,comm,
		     kcut,nloc,sz,data,pivot,cmp);

  if(debug) {
    integer x = 1;
    /* Barrier... */
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(x),&x,sum_integer);
    printf("{%s:%d} pid=%03d: parkstat returned, npivots = %d, pivot = %.3f\n",
	   __FILE__,__LINE__,
	   (int) pid,(int) npivots,((double *) pivot)[0]);
  }


  /* Local partitioning according to pivot */ {
    integer i = 0,j = nloc-1;
    
    while(i <= j) {
      if(cmp(idx_ptr(data,i,sz),pivot) <= 0) i++;
      else if(cmp(idx_ptr(data,j,sz),pivot) > 0) j--;
      else {
	swap(data,i,j,sz);
	i++;
	j--;
      }
    }
    nle = i;
    ngt = nloc - nle;
  }

  if(debug) printf("{%s:%d} pid=%03d: nloc = %d, nle = %d, ngt = %d, npivots = %d, pivot = %.3f\n",
	 __FILE__,__LINE__,(int) pid,(int) nloc,(int) nle,(int) ngt,(int) npivots,((double *) pivot)[0]);

  /*free(pivot);*/


  /* Count sizes of partitions (could be done by parkstat) */ {
    integer tmp[2];
    tmp[0] = nle;
    tmp[1] = ngt;
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   2,sizeof(integer),&tmp,sum_integer);
    nle_tot = tmp[0];
    ngt_tot = tmp[1];
  }

  if(debug) printf("{%s:%d} pid=%03d: nle_tot = %d, ngt_tot = %d\n",
	 __FILE__,__LINE__,(int) pid,(int) nle_tot,(int) ngt_tot);
  
  /* Balance keep and send partitions over processors, and exchange data */ {
    
    if(g_pid < pcut) {
      integer nle_left,ngt_left;

      /* Equalize large elements over left partition */
      ngt_left = equalize(pid0,pid0+pcut*pstride,pid,pstride,comm,
			  &ngt,nlocmax-nle,sz,idx_ptr(data,nle,sz));
      if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);

      if(debug) printf("{%s:%d} pid=%03d: (left) after equalize, ngt_left = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) ngt_left);
      if(debug) vprint(pid,nle,ngt,data);
      
      /* Send large elements to right partition */
      exchange(pid0,pid0+pcut*pstride,pid,pstride,
	       pid0+pcut*pstride,pid1,pstride,comm,
	       ngt_left,0,
	       ngt_tot-ngt_left,ngt_tot,
	       nlocmax-nle,sz,idx_ptr(data,nle,sz));

      /* Equalize small (remaining) elements over left partition */
      nle_left = equalize(pid0,pid0+pcut*pstride,pid,pstride,comm,
			 &nle,nlocmax,sz,data);

      if(debug) printf("{%s:%d} pid=%03d: (left) after equalize, nle_left = %d  nle = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) nle_left,(int) nle);
      if(debug) vprint(pid,nle,0,data);

      /* Receive small elemenets from right partition */
      exchange(pid0,pid0+pcut*pstride,pid,pstride,
	       pid0+pcut*pstride,pid1,pstride,comm,
	       nle_left,nle_tot,
	       nle_tot-nle_left,0,
	       nlocmax,sz,data);

      /* Deduce how many elements I have */ {
	integer a,b;
	partition(nle_tot,pcut,g_pid,&a,&b);
	nle = b-a;
      }

      if(debug) printf("{%s:%d} pid=%03d: (left) after recevie, nloc = nle = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) nle);


    } else {
      integer nle_right,ngt_right;

      if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);
      /* Equalize large elements over right partition */
      ngt_right = equalize(pid0+pcut*pstride,pid1,pid,pstride,comm,
			 &ngt,nlocmax-nle,sz,idx_ptr(data,nle,sz));

      if(debug) printf("{%s:%d} pid=%03d: (right) after equalize (1), ngt_right = %d  ngt = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) ngt_right,(int) ngt);

      /* Put small elements last */
      if(ngt <= nle)
	memswap(data,idx_ptr(data,nle,sz),ngt*sz);
      else
	memswap(data,idx_ptr(data,ngt,sz),nle*sz);

      
      /* Equalize small elements over right partition before receiving */
      nle_right = equalize(pid0+pcut*pstride,pid1,pid,pstride,comm,
	       &nle,nlocmax-ngt,sz,idx_ptr(data,ngt,sz));

      if(debug) printf("{%s:%d} pid=%03d: (right) after equalize (2), nle_right = %d  nle = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) nle_right,(int) nle);


      /* Put large elements back at the end */
      if(ngt <= nle)
	memswap(data,idx_ptr(data,nle,sz),ngt*sz);
      else
	memswap(data,idx_ptr(data,ngt,sz),nle*sz);

      /* Receive large elements from left partition */
      exchange(pid0+pcut*pstride,pid1,pid,pstride,
	       pid0,pid0+pcut*pstride,pstride,comm,
	       ngt_right,ngt_tot,
	       ngt_tot-ngt_right,0,
	       nlocmax-nle,sz,idx_ptr(data,nle,sz));
      /* Deduce how many large elements I have */ {
	integer a,b;
	partition(ngt_tot,np-pcut,g_pid-pcut,&a,&b);
	ngt = b-a;
      }

      if(debug) printf("{%s:%d} pid=%03d: (right) after receive, ngt = %d\n",
	     __FILE__,__LINE__,(int) pid,(int) ngt);


      /* Again, small elements at the end */
      if(nle <= ngt)
	memswap(data,idx_ptr(data,ngt,sz),sz*nle);
      else
	memswap(data,idx_ptr(data,nle,sz),sz*ngt);

      /* Send small elements to left partition */
      exchange(pid0+pcut*pstride,pid1,pid,pstride,
	       pid0,pid0+pcut*pstride,pstride,comm,
	       nle_right,0,
	       nle_tot-nle_right,nle_tot,
	       nlocmax-ngt,sz,idx_ptr(data,ngt,sz));

    }

  }

  if(serialize) tree_barrier(pid0,pid1,pid,pstride,comm);
  {
    integer ntot2 = (g_pid < pcut) ? nle : ngt;
    tree_allreduce(pid0,pid1,pid,pstride,comm,
		   1,sizeof(ntot2),&ntot2,sum_integer);
    if(g_pid == 0) if(ntot2 != ntot) {
	printf("Element count inconsistency ntot=%d ntot2=%d...",
	       (int) ntot,(int) ntot2);
	MPI_Abort(comm,5);
    }
  }

  STOPTIMER;
  if(g_pid < pcut) return nle;
  else return ngt;
}

integer balance(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc,integer nlocmax,integer sz,void *data,cmp_fun cmp) {
  STARTTIMER;
  void *pivot = malloc(sz);
  integer n = balance2(pid0,pid1,pid,pstride,comm,
		       nloc,nlocmax,sz,data,cmp,pivot);
  free(pivot);
  STOPTIMER;
  return n;
}
