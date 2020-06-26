/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef EXCHANGE__
#define EXCHANGE__

#include <mpi.h>
#include "integer.h"

#include "txmap.h"

void partition(integer n,integer np,integer pid,
	       integer idx0[1],integer idx1[1]);

void exchange_sched(integer pid0,integer pid1,integer pid,integer pstride,
		    integer rp0,integer rp1,integer rpstride,
		    MPI_Comm comm,
		    integer nhavetot,integer ngoaltot,
		    integer rnhavetot,integer rngoaltot,
		    txmap_struct *txmap);

void exchange(integer pid0,integer pid1,integer pid,integer pstride,
	      integer rp0,integer rp1,integer rpstride,
	      MPI_Comm comm,
	      integer nhavetot,integer ngoaltot,
	      integer rnhavetot,integer rngoaltot,
	      integer nlocmax,integer sz,void *data);

#endif

