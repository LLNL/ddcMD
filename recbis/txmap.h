/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef TXMAP__
#define TXMAP__

#include "integer.h"

typedef struct {
  integer pid,count;
} txmap_entry;

typedef struct {
  integer nused,nalloc;
  txmap_entry *list;
} txmap_struct;

#endif
