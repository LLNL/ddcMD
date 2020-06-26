/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef APPROX_MEDIAN__
#define APPROX_MEDIAN__

#include "treecomm.h"

#include "bisection_data.h"

double approx_kstat(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer nloc,bisection_data data[],int idir,double x0,double x1,
		    const double volume_cost,const double frac);
/*
Parallel approximate k-statistic computation, or sub-division according to cost function.
Elements are of type bisection_data (bisection_data.h), and represent 3D points. The cost
function is given by a weight for each point, plus a cost per unit length of the interval
[x0,x1]. The (approximate) sub-division point X is computed over x-, y- or z-coordinates.
The sub-division point is chosen so that `frac' of the total cost is to the left (< X),
and 1-`frac' to the right (> X). To perform k-statistic computation, set all point weights
to one, and set frac to k/ntot, where ntot is the totoal number f points across all
processors.

The problem is solved by a group of processors, with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , . . . , pid1).

Arguments:
   pid0    -- lowest processor MPI rank in group
   pid1    -- 1+higest processor MPI rank in group
   pid     -- MPI rank of this processor   
   pstride -- rank spacing between processor in group.
   comm    -- MPI communicator used for communication and which
               ranks pertain to

   nloc    -- number of elements on this processor
   data    -- pointer to elements
   idir    -- Cartesian direction in which to compute median, 0=x, 1=y, 2=z direction.
   x0      -- left-most (minimum) coordinate value in direction idir
   x0      -- right-most (maximum) coordinate value in direction idir
   volume_cost -- Cost per unit length in chosen direction
   frac        -- Return division point at this cost fraction

Return value:
    Returns an approximation of the point in `idir' direction X such that `frac' of the total
    cost is to the left (< X), and 1-`frac' to the right (> X). The accuracy is (x1-x0)/32^5.
    By changing niter and bins in approx_median.c, the accuracy can be improved or loosened.

The runtime is proportional to niter*(nloc + nbins*log(p), where p is the number of processors.
Defaults are niter = 5, and bins = 32.
*/





double approx_median(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc,bisection_data data[],int idir,double x0,double x1);
/*
Parallel approximate median computation. Elements are of type bisection_data (bisection_data.h),
and represent 3D points. The (approximate) median is computed over x-, y- or z-coordinates.


The problem is solved by a group of processors, with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , . . . , pid1).

Arguments:
   pid0    -- lowest processor MPI rank in group
   pid1    -- 1+higest processor MPI rank in group
   pid     -- MPI rank of this processor   
   pstride -- rank spacing between processor in group.
   comm    -- MPI communicator used for communication and which
               ranks pertain to

    nloc   -- number of elements on this processor
    data   -- pointer to elements
    idir   -- Cartesian direction in which to compute median, 0=x, 1=y, 2=z direction.
    x0     -- left-most (minimum) coordinate value in direction idir
    x0     -- right-most (maximum) coordinate value in direction idir

Return value:
    Returns an approximation of the median over 'idir'-coordinates. Accuracy is (x1-x0)/32^5.
    By changing niter and bins in approx_median.c, accuracy can be improved.

The runtime is proportional to niter*(nloc + nbins*log(p), where p is the number of processors.
Defaults are niter = 5, and bins = 32.
*/

#endif
