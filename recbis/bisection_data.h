/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef BISECTION_DATA__
#define BISECTION_DATA__

#include "integer.h"

typedef struct {
#ifdef BISECTION_COMPACT_SIZE
  float coord[3],weight;
  int origin,index_at_origin;
#else
  double coord[3],weight;
  integer origin,index_at_origin;
#endif
} bisection_data;

#endif
