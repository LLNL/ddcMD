/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef INTEGER__
#define INTEGER__

typedef long long int integer;

static inline void sum_integer(void *x,const void *dx) {
  ((integer *) x)[0] += ((const integer *) dx)[0];
}

static inline integer abs_integer(integer a) {
  if( a >= 0) return a;
  else return -a;
}

static inline integer min_integer(integer a,integer b) {
  if(a < b) return a;
  else return b;
}

static inline int cmp_integer(const void *ap,const void *bp) {
  integer
    a = ((const integer *) ap)[0],
    b = ((const integer *) bp)[0];
  if(a < b) return -1;
  else if(a == b) return 0;
  else return 1;
}

#endif
