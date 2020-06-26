#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "system.h"
#include "ddcenergy.h"

void forcetest(DDC *ddc, SYSTEM *sys)
{
	double *rx, *ry, *rz,x,y,z,delta=0.1250; 
	double *fx, *fy, *fz,nfx,nfy,nfz,ep,em; 
	double dfx,dfy,dfz; 
	int updateRate; 
	rx = sys->collection->state->rx; 
	ry = sys->collection->state->ry; 
	rz = sys->collection->state->rz; 
	fx = sys->collection->state->fx; 
	fy = sys->collection->state->fy; 
	fz = sys->collection->state->fz; 
	updateRate=ddc->updateRate; 
	ddc->updateRate=0; 
	ddcenergy(ddc, sys, 0);
   for (int k=0;k<15;k++) 
   {
   double err2Max = 0.0; 
   double err2 = 0.0; 
   int nsample = 100; 
   
	for (int j=0;j<nsample;j++)
	{
      int i = drand48()*sys->nlocal;
		x = rx[i];
		rx[i] = x + 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		ep =sys->energyInfo.eion; 
		rx[i] = x - 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		em =sys->energyInfo.eion; 
		rx[i]=x;
		nfx  =-(ep -em ) /delta; 

		y = ry[i];
		ry[i] = y + 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		ep =sys->energyInfo.eion; 
		ry[i] = y - 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		em =sys->energyInfo.eion; 
		ry[i]=y;
		nfy  =-(ep -em ) /delta; 

		z = rz[i];
		rz[i] = z + 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		ep =sys->energyInfo.eion; 
		rz[i] = z - 0.5*delta; 
		ddcenergy(ddc, sys, 0);
		em =sys->energyInfo.eion; 
		rz[i]=z;
		nfz  =-(ep -em ) /delta; 
		dfx = (nfx-fx[i]);
		dfy = (nfy-fy[i]);
		dfz = (nfz-fz[i]);
		double diff2 = dfx*dfx + dfy*dfy + dfz*dfz ;
      err2 += diff2; 
      if (diff2 > err2Max) err2Max = diff2; 
      /*
		if (dfx*dfx + dfy*dfy + dfz*dfz >1e-4) 
		{
			printf("%d %f %f %f ",i,nfx,nfy,nfz);
			printf("%f %f %f\n",fx[i],fy[i],fz[i]);
		}
      */
	}
	ddcenergy(ddc, sys, 0);
   double f2 =0; 
   for (unsigned i=0;i<sys->nlocal;i++) f2 += SQ(fx[i])+SQ(fy[i])+SQ(fz[i]); 
   double rmsF = sqrt(f2/nsample); 
   double rmsError=sqrt(err2/nsample);
   double errMax = sqrt(err2Max); 
   
   printf("Error: nsamples=%d delta=%.2e rmsF=%.2e rmsErr=%.2e  %.2e maxError=%.2e %.2e\n",nsample, delta, rmsF,rmsError,rmsError/rmsF,errMax,errMax/rmsF); 
  delta *= 0.5; 
   }
	ddc->updateRate=updateRate; 
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
