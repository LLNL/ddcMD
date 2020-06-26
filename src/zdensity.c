#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "analysis.h"

#include "three_algebra.h"
#include "object.h"

#include "simulate.h"
#include "system.h"
#include "ddcMalloc.h"
#include "external.h"
#include "io.h"
#include "box.h"
#include "units.h"
#include "mpiUtils.h"

enum SMEAR_METHOD
{
   SMEAR_IMPULSE, SMEAR_HAT
};


typedef struct zdensity_parms_st
{
   unsigned nz;
   char *filename; 
   double smearRadius;
   enum SMEAR_METHOD smearMethod;
}  ZDENSITY_PARMS ; 


ZDENSITY_PARMS* zdensity_parms(ANALYSIS *analysis)
{
   char* smearString; 
   
   ZDENSITY_PARMS* parms = analysis->parms= (ZDENSITY_PARMS *) ddcMalloc(sizeof(ZDENSITY_PARMS));
   OBJECT* obj = (OBJECT*) analysis;
   object_get(obj, "nz", &parms->nz, INT, 1, "0");
   object_get(obj, "filename", &parms->filename, STRING, 1, "zden.dat");
   object_get(obj, "smearRadius", &parms->smearRadius, WITH_UNITS, 1, "0", "l", NULL);
   object_get(obj, "smearMethod", &smearString, STRING, 1, "impulse");
   parms->smearMethod = SMEAR_IMPULSE;
   if (strcasecmp(smearString, "hat") == 0)
      parms->smearMethod = SMEAR_HAT;
   return parms; 
}

void zdensity_eval(ANALYSIS* analysis)
{
}

void zdensity_output(ANALYSIS *analysis)
{
   ZDENSITY_PARMS* parms = (ZDENSITY_PARMS*) analysis->parms; 
   SIMULATE* simulate = (SIMULATE*) analysis->parent; 
   SYSTEM* sys = simulate->system; 
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

   double density[parms->nz];
   for (unsigned ii=0; ii<parms->nz; ++ii)
      density[ii] = 0;
   
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
    wz[1]=0.0;
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
	 igz[1] = iWall.z;   if (igz[1] == (int)parms->nz) igz[1] = 0;
			
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
					
					
	 unsigned label = igz[kk];
	 if (label >= parms->nz)
	    label = parms->nz-1;
	 
	 assert(label < parms->nz);
	 density[label] += weight;
      }
		
   }

   // collect data from all tasks
   double totalDensity[parms->nz];
   MPI_Reduce(density, totalDensity, parms->nz, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);

   // write the data
   CreateSnapshotdir(simulate, NULL);
   if (getRank(0) == 0)
   {
      char filename[1024];
      snprintf(filename, 1024,"%s/%s", simulate->snapshotdir,parms->filename);
      FILE* file = fopen(filename, "w");
      double lc = units_convert(1.0, NULL, "Angstrom");
      double boxVol = box_get_volume(NULL)*lc*lc*lc;

      for (unsigned ii=0; ii<parms->nz; ++ii)
      {
	 double z = ( (ii+0.5) * (bbox.z/parms->nz) )/(bbox.z);
	 double dens = totalDensity[ii]*(parms->nz/boxVol);
	 fprintf(file, "%f %f %f\n", z, dens, totalDensity[ii]);
      }
      fclose(file);
   }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
