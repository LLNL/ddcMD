#include "loadBalance.h"

#include <string.h>
#include <assert.h>
#include <math.h>

#include "ddcMalloc.h"
#include "object.h"
#include "system.h"
#include "ddc.h"

// ToDo:
//
// 1.  Better default value for nz.  If left at zero the code will fail
// assertions when trying to findCenters.  This is a difficlt assertion
// failure to track down.


void domain_set(DDC* ddc);

enum SMEAR_METHOD
{
   SMEAR_IMPULSE, SMEAR_HAT
};

typedef struct zRampBalanceParms_st
{
   int nz;
   double smearRadius;
   enum SMEAR_METHOD smearMethod;
  
} ZRAMP_BALANCE_PARMS;

static void computeDensity(SYSTEM* sys, ZRAMP_BALANCE_PARMS* parms, double* density);
static void findCenters(double* density, unsigned nDensity, double* centers, unsigned nCenters);

void* zRampLoadBalance_init(LOAD_BALANCE* balancer)
{
   char* smearString; 
   
   ZRAMP_BALANCE_PARMS* parms = balancer->parms =
      (ZRAMP_BALANCE_PARMS*) ddcMalloc(sizeof(ZRAMP_BALANCE_PARMS));
   OBJECT* obj = (OBJECT*) balancer;
   object_get(obj, "nz", &parms->nz, INT, 1, "0");
   object_get(obj, "smearRadius", &parms->smearRadius, WITH_UNITS, 1, "0", "l", NULL);
   object_get(obj, "smearMethod", &smearString, STRING, 1, "impulse");
   parms->smearMethod = SMEAR_IMPULSE;
   if (strcasecmp(smearString, "hat") == 0)
      parms->smearMethod = SMEAR_HAT;
   
   return parms;
}

void zRampLoadBalanceFunction(LOAD_BALANCE* this, SYSTEM* sys, DDC* ddc)
{
   ZRAMP_BALANCE_PARMS* parms = this->parms; 

   double density[parms->nz];
   computeDensity(sys, parms, density);
   
   // computeDensity gives the particle density.  However, the work
   // density is proportional to the particle density squared.  We want
   // to use work density to find the domain centers.
   for (int ii=0; ii<parms->nz; ++ii)
      density[ii] *= density[ii];
   
   findCenters(density, parms->nz, ddc->dz, ddc->lz);
   
   domain_set(ddc);
}

void computeDensity(SYSTEM* sys, ZRAMP_BALANCE_PARMS* parms, double* density)
{
   STATE* state = sys->collection->state; 
   unsigned nlocal = sys->nlocal; 
   
   double* rz = state->rz; 
	
   THREE_VECTOR corner = box_get_corner(sys->box);
   THREE_VECTOR bbox = box_get_boundingbox(NULL);
   THREE_VECTOR deltai,scaled_corner;
   deltai.z = parms->nz/bbox.z; 
   scaled_corner.z = corner.z * deltai.z; 

   THREE_VECTOR lSmearInv = vzero;
   THREE_VECTOR lSmearHalf = vzero;
   if (parms->smearRadius > 0)
   {
      THREE_VECTOR lSmear;
      lSmear.z = MIN( 2.0*parms->smearRadius, bbox.z/(1.0*parms->nz));
      lSmearInv.z = 1.0/lSmear.z;
      lSmearHalf.z = 0.5 * lSmear.z;
   }

   double zBin[parms->nz];
   for (int ii=0; ii<parms->nz; ++ii) zBin[ii] = 0;

   for (unsigned iAtom=0; iAtom<nlocal; ++iAtom) 
   {
      THREE_VECTOR r;
      r.z=rz[iAtom]*deltai.z - scaled_corner.z; 

      unsigned  spread = 2;
      int igz[2];
      double wz[2];
      if ( parms->smearRadius <= 0)
      {
         spread = 1;
         igz[0] = r.z;
         igz[1] = 0.0; 
         wz[0]  = 1.0;
         wz[1]  = 0.0; 
      }
      else
      {
         spread = 2;
         THREE_INT iWall;
         iWall.z = floor(r.z + 0.5);

         THREE_VECTOR delta;
         delta.z = iWall.z - r.z;

         delta.z = MIN(delta.z, lSmearHalf.z);
         delta.z = MAX(delta.z, -lSmearHalf.z);

         igz[0] = iWall.z-1; if (igz[0] == -1)        igz[0] = parms->nz-1;
         igz[1] = iWall.z;   if (igz[1] == parms->nz) igz[1] = 0;

         switch (parms->smearMethod)
         {
            case SMEAR_IMPULSE:
               wz[0] = 0.5 + (delta.z * lSmearInv.z);
               break;
            case SMEAR_HAT:
               wz[0] = 0.5 + 2*delta.z*lSmearInv.z*(1.0 - fabs(delta.z)*lSmearInv.z);
               break;
            default:
               assert (1==0);
         }

         assert(wz[0] >= 0);
         assert(wz[0] <= 1);
         wz[1] = 1.0 - wz[0];
      }

      for (unsigned kk=0; kk<spread; ++kk)
      {
         double weight = wz[kk];
         if (weight < 1e-20)
            continue;

         int label = igz[kk];
         if (label >= parms->nz)
            label = parms->nz-1;

         assert(label < parms->nz);
         zBin[label] += weight;
      }

   }

   // collect data from all tasks
   MPI_Allreduce(zBin, density, parms->nz, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);

   double boxVol = box_get_volume(NULL);
   for (int ii=0; ii<parms->nz; ++ii)
   {
      //double z = ( (ii+0.5) * (bbox.z/parms->nz) )/(bbox.z);
      density[ii] *= (parms->nz/boxVol);
   }
}


void findCenters(double* density, unsigned nDensity, double* centers, unsigned nCenters)
{
   // integrate density
   double densitySum = 0;
   for (unsigned ii=0; ii<nDensity; ++ii)
      densitySum += density[ii];

   double densityTarget = densitySum/(1.0*nCenters);

   double walls[nCenters+1];
   walls[0] = 0;

   for (unsigned ii=0; ii<nCenters-1; ++ii)
   {
      double fPos = walls[ii];
      unsigned iPos = floor(walls[ii]);
      double sum = 0;
      double delta = 0;      
      while (1 == 1)
      {
         double dummy;
         double weight = 1.0 - modf(fPos, &dummy);

         if (sum + density[iPos]*weight > densityTarget)
            break;
         sum += density[iPos]*weight;
         delta += weight;

         ++iPos;
         fPos = iPos;
      }
      double fraction = (densityTarget - sum)/density[iPos];
      assert(fraction >= 0.0 && fraction <= 1.0);

      walls[ii+1] = walls[ii] + delta + fraction;
   }

   walls[nCenters] = nDensity;

   for (unsigned ii=0; ii<nCenters; ++ii)
   {
      int start = walls[ii];
      int end = walls[ii+1];
      double dummy;
      double weight = 1.0 - modf(walls[ii], &dummy);
      double sum = weight*density[start];
      for (int jj=start+1; jj<end; ++jj)
         sum += density[jj];
      weight = modf(walls[ii+1], &dummy);
      sum += weight*density[end];
      //      fprintf(ddcfile, "%u %f %f\n", ii, sum, densityTarget);
      //      fflush(ddcfile);
      assert(fabs(sum - densityTarget) < 0.01);
   }

   centers[1] = 0.5*(walls[1]+walls[2]);
   centers[0] = 2.0*walls[1] - centers[1];
   for (unsigned ii=2; ii<nCenters; ++ii)
   {
      centers[ii] = walls[ii] + (walls[ii] - centers[ii-1]);
      assert(centers[ii] < walls[ii+1]);
   }
   for (unsigned ii=0; ii<nCenters; ++ii)
      centers[ii] /= (nDensity*1.0);
}

