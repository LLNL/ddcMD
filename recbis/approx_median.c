/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <stdlib.h>
#include "treecomm.h"
#include "approx_median.h"

#include <stdio.h>
#include <math.h>

#include "rectimer.h"
#include "functimer.h"

const int nbins = 32,niter = 5;

void sum_hist(void *x,const void *dx) {
  double *hist = (double *) x;
  const double *dhist = (const double *) dx;
  int i;
  for(i = 0; i<nbins; i++)
    hist[i] += dhist[i];
}

double approx_median(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer nloc,bisection_data data[],int idir,double x0,double x1) {
  return approx_kstat(pid0,pid1,pid,pstride,comm,
		      nloc,data,idir,x0,x1,
		      /* volume_cost: */ 0.0,
		      /* frac:        */ 0.5);
}

double approx_kstat(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer nloc,bisection_data data[],int idir,double x0,double x1,
		    const double volume_cost,const double frac) {
  STARTTIMER;
  const int dbg = 0;
  const double range0 = x1-x0;
  double wleft = 0.0,wright = 0.0;

  int ptrs[nbins],*next = (int *) malloc(nloc * sizeof(*next));
  int lastbin = -1;
  int iter;

  for(iter = 0; iter<niter; iter++) {
    double hist[nbins],wtot;

    {
      int j;
      for(j = 0; j<nbins; j++) hist[j] = 0.0;
    }

    if(lastbin < 0) {
      /* First iteration: Put particles in respective bins */
      if(dbg && pid == 0)
	printf("Binning, first iteration x0=%.3f, x1=%.3f, frac=%.3f\n",x0,x1,frac);
      {
	int j;
	for(j = 0; j<nbins; j++) ptrs[j] = -1;
      }
      {
	int j;
	for(j = 0; j<nloc; j++) {
	  int ibin = (data[j].coord[idir] - x0)/(x1-x0) * nbins;
	  if(ibin < 0) ibin = 0;
	  else if(ibin >= nbins) ibin = nbins - 1;
	  hist[ibin] += data[j].weight;
	  next[j] = ptrs[ibin];
	  ptrs[ibin] = j;
	}
      }
    } else {
      /* Distribute particles in lastbin to new binset */
      int p = ptrs[lastbin];
      if(dbg && pid == 0)
	printf("Binning, subsequent iteration (iter = %d)\n",iter);
      {
	int j;
	for(j = 0; j<nbins; j++) ptrs[j] = -1;
      }
      while(p >= 0) {
	int ibin = (data[p].coord[idir] - x0)/(x1-x0) * nbins;
	if(ibin < 0) ibin = 0;
	else if(ibin >= nbins) ibin = nbins - 1;
	hist[ibin] += data[p].weight;
	
	/* Move particle from lastbin to current bin */ {
	  const int np = next[p];
	  next[p] = ptrs[ibin];
	  ptrs[ibin] = p;
	  p = np;
	}
      }
    }

    if(dbg && pid == 0) printf("Entering tree_allreduce\n");

    /* Redeuce histogram over all tasks */ {
      //MPI_Allreduce(hist,nbins,MPI_DOUBLE,MPI_SUM,comm);
      tree_allreduce(pid0,pid1,pid,1,comm,
		     1,sizeof(hist),hist,sum_hist);

      /* Increment histogram with volume_cost per bin */ {
	const double bin_cost = volume_cost*(x1-x0)/(range0*nbins);
	{
	  int j;
	  for(j = 0; j<nbins; j++) hist[j] += bin_cost;
	}
      }

    }

    /* Compute total weight, could be done once only */ {
      wtot = wleft + wright;
      {
	int j;
	for(j = 0; j<nbins; j++)
	  wtot += hist[j];
      }
    }
    if(dbg && pid == 0) printf("wtot = %.3f, wleft = %.3f, wright = %.3f\n",
			       wtot,wleft,wright);

    if(dbg && pid == 0) printf("Find median bin\n");
    double wsum = wleft;
    /* Find median bin, i.e. first bin where cumulative weight is more than half the total */ {
      int j = 0;
      double /*wsum = 0.0,*/nx0,nx1;

      double wsummax = wsum;
      {
	int jj;
	for(jj = 0; jj<nbins; jj++) wsummax += hist[j];
      }
      if(wsummax-wsum < 1e-15*nbins*wsummax) j = (nbins+1)/2;
      else {
	while(wsum < frac*wtot) {
	  if(j >= nbins) {
	    //printf("Warning: wsum[nbins-1] < frac*wtot\n");
	    j = nbins-1;
	    break;
	  }
	  
	  wsum += hist[j++];
	  if(dbg && pid == 0)
	    printf("j = %d, wsum = %.3f, wtot = %.3f, frac*wtot = %.3f\n",
		   j,wsum,wtot,frac*wtot);
	}
      }
      {
	int i;
	for(i = 0; i<j-1; i++) wleft += hist[i];
	for(i = j; i<nbins; i++) wright += hist[i];
      }
      nx0 = x0 + (j-1)*(x1-x0)/nbins;
      nx1 = x0 + j*(x1-x0)/nbins;
      x0 = nx0;
      x1 = nx1;
      lastbin = j-1;
    }

    if(dbg && pid == 0)
      printf("Median bin is lastbin=%d, x0=%.5f, x1 = %.5f (wsum = %.3f, wtot = %.3f)\n",
	     lastbin,x0,x1,wsum,wtot);
  }

  free(next);

  STOPTIMER;
  /* Define median as halfway between bin boundaries */ {
    const double xcut = 0.5*(x1 + x0);
    return xcut;
  }
}

#ifdef BUILD_MAIN
#include <stdlib.h>
#include <stdio.h>

#include "back_communicate.h"

int bisec_cmp(const void *a,const void *b) {
  const double
    xa = ((const bisection_data *) a)->coord[0],
    xb = ((const bisection_data *) b)->coord[0];
  if(xa < xb) return -1;
  else if(xb < xa) return 1;
  else return 0;
}
integer get_dest(const void *a) {
  return ((const bisection_data *) a)->dest;
}

int main(int argc,char *argv[]) {
  const int dbg = 0;
  int nalloc,nmax = atoi(argv[1]),niter = atoi(argv[2]);
  bisection_data *data;
  int pid,np;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&pid);
  MPI_Comm_size(comm,&np);
  
  //srand((331 * pid) ^ 4723681);

  if(pid == 0) {
    printf("Using nmax = %d, niter = %d\n",
	   nmax,niter);
    fflush(stdout);
  }

  nalloc = 2*nmax;
  data = (bisection_data *) malloc(nalloc * sizeof(bisection_data));
  /* Check memory allocation error... */ {
    int flag = (data == NULL),errcnt;
    MPI_Allreduce(&flag,&errcnt,1,MPI_INT,MPI_SUM,comm);
    if(errcnt > 0) {
      if(pid == 0) {
	printf("Error allocating memory for %d elements: %d tasks failed.\n",
	       nalloc,errcnt);
	fflush(stdout);
      }
      MPI_Finalize();
      exit(1);
    }
  }

  for(int iter = 0; iter<niter; iter++) {
    int n;
    double xcut,emed;

    if(pid == 0) {
      n = rand() % nmax;

      if(dbg) {
	printf("Iteration %d, total number of elements is %d, task average is %.3f\n",
	       iter+1,n,n/(double) np);
	fflush(stdout);
      }

      for(int i = 0; i<n; i++) {
	data[i].dest = rand() % np;
	for(int j = 0; j<3; j++)
	  data[i].coord[j] = rand() / (double) RAND_MAX;
	data[i].weight = rand() / (double) RAND_MAX;
      }

      if(dbg) {
	printf("Sorting...\n");
	fflush(stdout);
      }
      qsort(data,n,sizeof(*data),bisec_cmp);

      /* Compute weighted median */ {
	double wsum = 0.0,wtot = 0.0;
	int j = 0;
	for(int i = 0; i<n; i++)
	  wtot += data[i].weight;
	while(wsum < 0.5*wtot) wsum += data[j++].weight;
	emed = data[j-1].coord[0];
	if(dbg) {
	  printf("Exact median is between %15.10f and %15.10f\n",
		 data[j-1].coord[0],data[j].coord[0]);
	  fflush(stdout);
	}
      }
    } else n = 0;
    
    MPI_Barrier(comm);
    if(dbg && pid == 0) {
      printf("Doing alltoall_sparse()\n");
      fflush(stdout);
    }
    MPI_Barrier(comm);

    n = alltoall_sparse(0,np,pid,1,comm,n,nalloc,sizeof(*data),data,get_dest);

    MPI_Barrier(comm);
    if(dbg && pid == 0) {
      printf("Doing approx_median()\n");
      fflush(stdout);
    }
    MPI_Barrier(comm);


    xcut = approx_median(0,np,pid,1,comm,
			 n,data,0,0.0,1.0);
    MPI_Bcast(&emed,1,MPI_DOUBLE,0,comm);
    if(fabs(xcut - emed) > 1e-14*fabs(emed)) {
      if(pid == 0) {
	printf("Iteration = %d, n = %d\n",iter,n);
	printf("xcut = %30.20e\nemed = %30.20e\ndiff = %30.20e\n",
	       xcut,emed,xcut-emed);
	fflush(stdout);
      }
      MPI_Finalize();
      exit(1);
    }

    MPI_Barrier(comm);
    if(dbg && pid == 0) {
      printf("Approximate median is %15.10f\n",xcut);
      fflush(stdout);
    }
    MPI_Barrier(comm);
  }

  free(data);
  MPI_Finalize();
  return 0;
}

#endif
