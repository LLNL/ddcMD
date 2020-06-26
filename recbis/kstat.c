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

#include "kstat.h"
#include "genptr.h"
#include "integer.h"

static void vprint(void *v_in,int n) {
  int *v = (int *) v_in;
  int i;
  for(i = 0; i<n; i++)
    printf("  %10d",v[i]);
  printf("\n");
}

void kstat(integer k,integer n,integer sz,void *v,cmp_fun cmp,void *pivot) {
  integer i;

  const int debug = 0;

  for(i = 0; i<n; i+=5) {
    integer m = min_integer(n-i,5);
    if(debug)
      qsort(idx_ptr(v,i,sz),m,sz,cmp);
    else {
      integer j;
      for(j = 1; j<m; j++) {
	integer k = j;
	while(k > 0 && cmp(idx_ptr(v,i+k,sz),idx_ptr(v,i+k-1,sz)) < 0) {
	  swap(v,i+k,i+k-1,sz);
	  k--;
	}
      }

      for(j = 1; j<m; j++)
	if(cmp(idx_ptr(v,i+j-1,sz),idx_ptr(v,i+j,sz)) > 0) {
	  printf("Sorting error. j=%lld, m=%lld\n",
		 (long long int) j,(long long int) m);
	  //vprint(v+i,m);
	  exit(1);
	}
    }
  }

  if(n <= 5)
    memcpy(pivot,idx_ptr(v,k,sz),sz); /* kval = v[k]; */
  else {
    for(i = 0; i<n/5; i++) {
      integer idx = 5*i + 2;
      swap(v,idx,i,sz);
    }

    {
      integer m = n % 5;
      if(m > 0) {
	integer idx = n - m/2 - 1;
	swap(v,idx,i,sz);
	i++;
      }
    }
  
    {
      integer j;
      void *kval_check = NULL;

      integer npiv,nlt,ngt;
      integer npiv_c,nlt_c,ngt_c;
      integer i0,j1;

      kstat((i-1)/2,i,sz,v,cmp,pivot);

      if(debug) {
	void *v2 = (void *) malloc(sz * n);
	memcpy(v2,v,n*sz);
	qsort(v2,n,sz,cmp);

	kval_check = (void *) malloc(sz);
	memcpy(kval_check,idx_ptr(v2,k,sz),sz);

	i = 0;
	while(i < n && cmp(idx_ptr(v2,i,sz),pivot) < 0) i++;
	j = i;
	while(j < n && cmp(idx_ptr(v2,j,sz),pivot) == 0) j++;

	npiv_c = j-i;
	nlt_c = i;
	ngt_c = n-j;

	free(v2);
      }

      i = 0; i0 = 0;
      j = n-1; j1 = n;

      while(i <= j) {
	int cond_i,cond_j;
	if((cond_i = cmp(idx_ptr(v,i,sz),pivot)) == 0 /* v[i] == p */) {
	  if(i > i0) {
	    swap(v,i,i0,sz);
	  }
	  i0++;
	  i++;
	} else if((cond_j = cmp(idx_ptr(v,j,sz),pivot)) == 0 /* v[j] == p */) {
	  j1--;
	  if(j < j1) {
	    swap(v,j,j1,sz);
	  }
	  j--;
	} else if(cond_i < 0 /* v[i] < p */)
	  i++;
	else if(cond_j > 0 /* v[j] > p */)
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

      if(debug) {
	if(npiv != npiv_c ||
	   nlt != nlt_c  || 
	   ngt != ngt_c) {
	  printf("Partition error... "
		 "npiv=%lld,npiv_c=%lld,nlt=%lld,nlt_c=%lld,ngt=%lld,ngt_c=%lld,"
		 "n=%lld,p=%lld,k=%lld\n",
		 (long long int) npiv,(long long int) npiv_c,
		 (long long int) nlt,(long long int) nlt_c,
		 (long long int) ngt,(long long int) ngt_c,
		 (long long int) n,(long long int) *(int *) pivot,(long long int) k);
	  vprint(v,n);
	  exit(1);
	}
      }


      if(k < nlt)
	kstat(k,nlt,sz,idx_ptr(v,i0,sz),cmp,pivot);
      else if(k < nlt+npiv) {
	/* kval = p; */
	/* pivot already set... */

	/*
	printf("pivot already set (n=%d,k=%d,pivot=%d\n",
	       (int) n,(int) k,*(int *) pivot);
	*/
      } else
	kstat(k-nlt-npiv,ngt,sz,idx_ptr(v,i,sz),cmp,pivot);

      if(debug) {
	if(cmp(pivot,kval_check) != 0 /* kval != kval_check */) {

	  printf("Broken... kval=%d, kval_check=%d\n",
		 *(int *) pivot,*(int *) kval_check);

	  printf("Broken...\n");
	  exit(1);
	}
	free(kval_check);
      }
    }
  }
  /* return kval; */
}
