/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef GENPTR__
#define GENPTR__

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <alloca.h>

#include "integer.h"

typedef int (*cmp_fun)(const void *a,const void *b);

static inline void *idx_ptr(void *v,integer idx,integer sz) {
  long long int p = (long long int) v;
  assert(sizeof(p) >= sizeof(v));
  return (void *) (p + idx*sz);
}

static inline void swap(void *v,integer i,integer j,integer sz) {
  const int max_stack_buf = 512;
  int ibuf[8];
  void *ptr;

  if(sz <= (integer) sizeof(ibuf)) ptr = (void *) ibuf;
  else if(sz <= max_stack_buf) ptr = alloca(sz);
  else ptr = malloc(sz);

  memcpy(ptr,idx_ptr(v,i,sz),sz);
  memcpy(idx_ptr(v,i,sz),idx_ptr(v,j,sz),sz);
  memcpy(idx_ptr(v,j,sz),ptr,sz);

  if(sz > max_stack_buf) free(ptr);
}

static inline void memswap(void *a,void *b,int sz) {
  const int bufsize = 512;
  void *t = alloca(bufsize);
  int i;

  for(i = 0; i<sz; i+=bufsize) {
    int n = sz-i;
    if(n > bufsize) n = bufsize;
    memcpy(t,a,n);
    memcpy(a,b,n);
    memcpy(b,t,n);
    a = idx_ptr(a,1,n);
    b = idx_ptr(b,1,n);
  }
}

#endif

