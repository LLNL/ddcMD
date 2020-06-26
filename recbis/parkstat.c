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

#include "parkstat.h"

#include "kstat.h"
#include "treecomm.h"

#include "genptr.h"
#include "integer.h"

/* Only used by qsort for debugging.. */
static int cmp_int(const void *ap,const void *bp) {
  int
    a = ((const int *) ap)[0],
    b = ((const int *) bp)[0];
  if(a < b) return -1;
  else if(a == b) return 0;
  else return 1;
}

integer parkstat(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
	      integer k,integer nloc,integer sz,void *vloc,void *pivot,cmp_fun cmp) {

  const int debug = 0;
  const integer gsize_max = 20;
  const integer np = (pid1-pid0)/pstride;
  const integer gid = ((pid-pid0)/pstride) / gsize_max;
  const integer g_pid = ((pid-pid0)/pstride) % gsize_max;
  const integer gsize = min_integer(np - gsize_max*gid,gsize_max);
  const integer g_pid0 = pid - pstride*g_pid;
  const integer g_pid1 = g_pid0 + pstride*gsize;
  const integer ngroups = (np+gsize_max-1)/gsize_max;

  integer nget,kval;

  integer junk = 1;

  void *vrem = (void *) malloc(sz * gsize_max);


  if(debug) {
    integer ntot = nloc;
    tree_allreduce(pid0,pid1,pid,pstride,comm,1,sizeof(ntot),&ntot,sum_integer);
    if(pid == 0) printf("  Gather phase. k=%d, nloc=%d, ntot=%d\n",
			(int) k,(int) nloc,(int) ntot);
    if(k >= ntot) {
      printf("pid=%d Error, dying, k=%d > ntot=%d\n",(int) pid,(int) k,(int) ntot);
      exit(1);
    }
  }

  if(nloc > 0) {
    if(debug) {
      int x;
      int *y = (int *) malloc(sizeof(int) * nloc);

      printf("Calling kstat... nloc=%d k=%d\n",(int) nloc,(int) ((nloc-1)/2));

      for(x = 0; x<nloc; x++)
	y[x] = ((int *) vloc)[x];
      qsort(y,nloc,sizeof(int),cmp_int);
      for(x = 0; x<nloc; x++)
	printf("%4d%s ",(int) y[x],(x==(nloc-1)/2) ? "*":" ");
      printf("\n");
      free(y);
    }
    kstat((nloc-1)/2,nloc,sz,vloc,cmp,vrem);

    if(debug)
      printf("kstat returned vrem[0]=%d\n",((int *) vrem)[0]);
    nget = tree_gather_inc(g_pid0,g_pid1,pid,pstride,comm,
			   1,1,sz,vrem);
    if(debug)
      printf("tree_grather_inc(1) returned nget=%d\n",(int) nget);
  } else {
    nget = tree_gather_inc(g_pid0,g_pid1,pid,pstride,comm,
			   1,0,sz,vrem);
    if(debug)
      printf("tree_grather_inc(0) returned nget=%d\n",(int) nget);
  }

  if(debug) {
    tree_allreduce(pid0,pid1,pid,pstride,comm,1,sizeof(junk),&junk,sum_integer);
    if(pid == 0) printf("  Recursive reduction + broadcast. k=%d, nloc=%d\n",(int) k,(int) nloc);
  }

  if(g_pid == 0) {
    const integer g_stride = pstride*gsize_max;
    integer ngel = (nget > 0) ? 1:0;
    void *ploc;
    tree_allreduce(pid0,pid0+ngroups*g_stride,pid,g_stride,comm,
		   1,sizeof(ngel),&ngel,sum_integer);

    ploc = malloc(sz);
    if(nget > 0)
      kstat((nget-1)/2,nget,sz,vrem,cmp,ploc);
    if(ngroups > 1) {
      parkstat(pid0,pid0+ngroups*g_stride,pid,g_stride,comm,
	       (ngel-1)/2,(nget > 0) ? 1:0,sz,ploc,pivot,cmp);
    } else {
      memcpy(pivot,ploc,sz);
    }
    free(ploc);
  }
  free(vrem);
  tree_broadcast(g_pid0,g_pid1,pid,pstride,comm,
		 1,sz,pivot);

  if(debug) {
    tree_allreduce(pid0,pid1,pid,pstride,comm,1,sizeof(junk),&junk,sum_integer);
    if(pid == 0) printf("  Partition phase. pivot=%d, nloc=%d\n",((int *) pivot)[0],(int) nloc);
  }

  /* Partition local data, get global partition size, and recurse */ {
    const void *p = pivot;
    const integer n = nloc;
    void *v = vloc;
    integer npiv,nlt,ngt;
    integer i,j,i0,j1;

    i = 0; i0 = 0;
    j = n-1; j1 = n;
    
    while(i <= j) {
      int ci,cj;
      if( (ci = cmp(idx_ptr(v,i,sz),p)) == 0 /*v[i] == p*/) {
	if(i > i0) swap(v,i,i0,sz);
	i0++;
	i++;
      } else if( (cj = cmp(idx_ptr(v,j,sz),p)) == 0 /*v[j] == p*/) {
	j1--;
	if(j < j1) swap(v,j,j1,sz);
	j--;
      } else if(ci < 0)
	i++;
      else if(cj > 0 /*v[j] > p*/)
	j--;
      else {
	swap(v,i,j,sz);
	i++;
	j--;
      }
    }
    
    npiv = i0 + (n-j1);
    nlt = i-i0;
    ngt = j1-j-1;

    /* Find global partition sizes, and recurse */ {
      integer npiv_tot,nlt_tot,ngt_tot,n_tot;
      integer arr[4];
      
      arr[0] = npiv;
      arr[1] = nlt;
      arr[2] = ngt;
      arr[3] = nloc;
      tree_allreduce(pid0,pid1,pid,pstride,comm,
		     sizeof(arr)/sizeof(integer),sizeof(integer),arr,sum_integer);
      npiv_tot = arr[0];
      nlt_tot = arr[1];
      ngt_tot = arr[2];
      n_tot = arr[3];

      if(debug) {
	integer kk;
	for(kk = pid0; kk<pid1; kk+=pstride) {
	  tree_allreduce(pid0,pid1,pid,pstride,comm,
			 1,sizeof(junk),&junk,sum_integer);

	  if(pid == kk)
	    printf("pid=%d  Recursion phase. k=%d, nlt_tot=%d npiv_tot=%d ngt_tot=%d n_tot=%d, "
		   "nlt=%d npiv=%d ngt=%d nloc=%d\n",
		   (int) pid,(int) k,
		   (int) nlt_tot,(int) npiv_tot,(int) ngt_tot,(int) n_tot,
		   (int) nlt,(int) npiv,(int) ngt,(int) nloc);
	  //usleep(100000);
	}
      }
      if(k < nlt_tot)
	kval = parkstat(pid0,pid1,pid,pstride,comm,
			k,nlt,sz,idx_ptr(v,i0,sz),pivot,cmp);
      else if(k < nlt_tot+npiv_tot) {
	kval = npiv_tot;
      } else
	kval = parkstat(pid0,pid1,pid,pstride,comm,
			k-nlt_tot-npiv_tot,ngt,sz,idx_ptr(v,i,sz),pivot,cmp);
    }
  }

  if(debug) {
    tree_allreduce(pid0,pid1,pid,pstride,comm,1,sizeof(junk),&junk,sum_integer);
    if(pid == 0) printf("Returning with kval=%d\n",(int) kval);
  }

  return kval;
}
