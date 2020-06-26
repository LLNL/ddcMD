/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef KSTAT__
#define KSTAT__

#include "integer.h"
#include "genptr.h"

void kstat(integer k,integer n,integer sz,void *v,cmp_fun cmp,void *pivot);

/*
Compute the k-statistic of a data set, meaning if the data set
were to be sorted in increasing order, return the k'th element.

This routine runs in time proportional to n (linear time).

Arguments:
   k     -- which element to return
   n     -- number of elements in input array
   sz    -- size of each element (e.g. as from sizeof() operator)
   v     -- pointer to elements
   cmp   -- function pointer to comparison function
   pivot -- pointer to element, which on return will be set to k-statsitic.


This routine runs in time proportional to n (linear time). It is assumed
that cmp takes constant time.
*/

#endif
