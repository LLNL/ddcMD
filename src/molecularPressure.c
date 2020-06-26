#include "three_algebra.h"
#include "system.h"
#include "collection.h"
#include "state.h"
#include "preduce.h"
#include "molecule.h"
#include "codata.h"
THREE_SMATRIX molecularVirialOld(SYSTEM *sys,THREE_SMATRIX virial)
{
   STATE *state = sys->collection->state; 
   THREE_VECTOR delta[sys->nlocal];   //Thought needed for parallelization. 
   for(unsigned i=0;i<sys->nlocal;i++) delta[i]=vzero; 
   moleculeDelta(sys->moleculeClass, state, delta);
   for (unsigned i=0;i<sys->nlocal;i++)
   {
      virial.xx -=  delta[i].x*state->fx[i]; 
      virial.yy -=  delta[i].y*state->fy[i]; 
      virial.zz -=  delta[i].z*state->fz[i]; 
   }
   return virial; 
}
THREE_SMATRIX  molecularVirial(SYSTEM *sys, THREE_SMATRIX virial)
{
   MOLECULE *molecule =  sys->moleculeClass->molecule; 
   STATE *state = sys->collection->state; 
   int nMolecules =  sys->moleculeClass->nMolecules; 
   for (int k=0;k<nMolecules;k++)
   {
      MOLECULETYPE  *type = molecule[k].type; 
      int a = type->ownershipSpeciesOffset; 
      int i0 = molecule[k].list[a];
      THREE_VECTOR r0 = {state->rx[i0],state->ry[i0],state->rz[i0]};
      double M=0.0; 
      THREE_VECTOR R=vzero; 
      THREE_VECTOR d[type->nSpecies]; 
      for (int a=0;a<type->nSpecies;a++)
      {
         int i = molecule[k].list[a];
         double mass = ((ATOMTYPE_PARMS *)(state->species[i]->parm))->mass;
         VSET(d[a],state->rx[i]-r0.x,state->ry[i]-r0.y,state->rz[i]-r0.z); 
         nearestImage(&d[a].x,&d[a].y,&d[a].z); 
         VSVOP(R,+=,mass,*,d[a]);
         M += mass; 
      }
      VSCALE(R,1.0/M); 
      for (int a=0;a<type->nSpecies;a++)
      {
         VOP1(d[a],-=,R);
         int i = molecule[k].list[a];
         virial.xx -=  d[a].x*state->fx[i]; 
         virial.yy -=  d[a].y*state->fy[i]; 
         virial.zz -=  d[a].z*state->fz[i]; 
      }
   }
   return virial; 
}
THREE_SMATRIX  molecularPressure(SYSTEM *sys, THREE_SMATRIX virial, double T)
{
   int N = sys->moleculeClass->nMoleculesGlobal; 
   double vol = box_get_volume(NULL); 
   THREE_SMATRIX  pTensor = molecularVirial(sys,virial);
   pTensor.xx += N*kB*T;
   pTensor.yy += N*kB*T;
   pTensor.zz += N*kB*T;
   SMATNORM(pTensor, vol);
   return pTensor; 
}
