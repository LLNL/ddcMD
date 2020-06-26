#include "langevin.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "system.h"
#include "ddcMalloc.h"
#include "random.h"
#include "eq.h"
#include "codata.h"
#include "units.h"
#include "langevin.h"
#include "hycop.h"

void skupskyParse(GROUP *group);
void localYukawaParse(GROUP *group);
void hycopLangevinParse(GROUP *group);


void langevin_write_dynamics(GROUP*g, FILE *file)
{
   LANGEVIN_PARMS *p; 
   p=g->parm ;
   fprintf(file,"%s %s { Teq=%f ;}\n",g->name,g->objclass,units_convert(p->Teq,NULL,"T"));
}
double langevin_getTemperature(LANGEVIN_PARMS *parms, double time)
{
   SYSTEM *system = system_getSystem(NULL); 
   double Teq=0; 
   if (parms->Teq_dynamics == GLOBAL_ENERGY ) 
   {
      double energy; 
      system_get(system,ENERGY,(void *)&energy);
      gid_type  nglobal = system_getNglobal(system);
      if (!parms->total_energy_set) 
      {
         parms->total_energy=parms->Teq*parms->Cp*nglobal+energy; 
         parms->total_energy_set = 1; 
      }
      Teq=(parms->total_energy - energy )/(parms->Cp * nglobal); 
   }
   if (parms->Teq_dynamics == EXPLICIT_TIME )
   {
      Teq = (parms->Teq_vs_time->function) (time, (void *)parms->Teq_vs_time);
   }
   return Teq;
}

void langevin_Update(GROUP *group, int mode, STATE *state, double time, double dt)
{
   LANGEVIN_PARMS *parms=group->parm ;
   parms->Teq = langevin_getTemperature(parms,time); 
   parms->kBT[0] = kB*parms->Teq; 
   return; 
}
static void *normalParse(GROUP *group)
{
   OBJECT *obj = (OBJECT*)group; 
   LANGEVIN_PARMS *parms = group->parm; 
   parms->type = NORMALLANGEVIN; 
   parms->strideTau =0; 
   parms->stridekBT =0; 
   parms->tau = malloc(sizeof(double)); 
   parms->kBT = malloc(sizeof(double)); 
   char *string; 
   object_get(obj, "Teq_dynamics", &string, STRING, 1, "EXPLICIT_TIME");
   if (!strcmp(string,"EXPLICIT_TIME"))  parms->Teq_dynamics = EXPLICIT_TIME;
   if (!strcmp(string,"GLOBAL_ENERGY"))  parms->Teq_dynamics = GLOBAL_ENERGY;
   if (parms->Teq_dynamics == GLOBAL_ENERGY) 
   {
      object_get(obj, "Cp", &parms->Cp, WITH_UNITS, 1, "1.0","m*l^2/t^2/T",NULL);
      object_get(obj, "Teq", &parms->Teq, WITH_UNITS, 1, "0.0","T",NULL);
      parms->total_energy_set=0; 
   }
   if (parms->Teq_dynamics == EXPLICIT_TIME) 
   {
      group->write_dynamics = langevin_write_dynamics;
      object_get(obj, "Teq", &string, LITERAL, 1, NULL);
      parms->Teq_vs_time = eq_parse(string,"T","t");
      parms->Teq = (double )((parms->Teq_vs_time->function) (0.0, (void *)parms->Teq_vs_time));
   }
   free(string); 
   object_get(obj, "tau", parms->tau, WITH_UNITS, 1, "1.0","t",NULL);

   return parms; 
}
void langevin_velocityUpdate(int mode, int  k,GROUP *group, STATE *state, double time, double dt)
{  
   LANGEVIN_PARMS *p=group->parm ; 
   RANDOM *random=p->random; 
   void *randomParms = random_getParms(random, k); 
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 
   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;

   THREE_VECTOR   v = p->vcm; 

   double kBT = p->kBT[k*p->stridekBT]; 
   double tau = p->tau[k*p->strideTau]; 
   double a = exp(-dt/tau);
   double c = dt/mass; 
   double d = sqrt(2.0*dt*kBT/(mass*tau));
   THREE_VECTOR g;
   g = gasdev3d(random,randomParms); 
   switch(mode)
   {
      case FRONT_TIMESTEP:
         vx[k]= v.x + a*(vx[k]-v.x) + c*fx[k] + d*g.x;
         vy[k]= v.y + a*(vy[k]-v.y) + c*fy[k] + d*g.y;
         vz[k]= v.z + a*(vz[k]-v.z) + c*fz[k] + d*g.z;
         break;
      case BACK_TIMESTEP:
         vx[k]= v.x + a*((vx[k]-v.x) + c*fx[k] + d*g.x);
         vy[k]= v.y + a*((vy[k]-v.y) + c*fy[k] + d*g.y);
         vz[k]= v.z + a*((vz[k]-v.z) + c*fz[k] + d*g.z);
         break;
   }
}
void langevin_velocityUpdateKernel(int mode, GROUP *group, STATE *state, double time, double dt)
{
   for (int i=0;i<state->nlocal;i++)
   {
      if (state->group[i] != group) continue; 
      langevin_velocityUpdate(mode,i,group,state,time,dt); 
   }
}
void langevin_parms(GROUP *group)
{
   group->itype = LANGEVIN;
   group->write_dynamics = NULL;
   group->velocityUpdate= (void (*)(int, int, GROUP *, void *,double, double))langevin_velocityUpdate; 
   group->Update= (void (*) (GROUP*, int, void *, double, double))langevin_Update; 
   group->field_names=""; 

   OBJECT *obj = (OBJECT*)group; 


   char *string,*end; 
   object_get(obj, "tau", &string, STRING, 1, "NotDefined");
   strtod(string,&end);  
   if (string != end) {free(string); string = strdup("NORMAL");}   //check to see if tau is a numeric value; 
   if (!strcasecmp(string,"SKUPSKYLANGEVIN"))  skupskyParse(group);
   if (!strcasecmp(string,"LOCALYUKAWA"))  localYukawaParse(group);
   if (!strcasecmp(string,"HYCOP"))  hycopLangevinParse(group);
  //HYCOPLANGEVIN_PARMS *parmHycop = (HYCOPLANGEVIN_PARMS *)group->parm;
   if (!strcasecmp(string,"NORMAL"))  
   {
      group->parm = ddcMalloc(sizeof(LANGEVIN_PARMS));
      normalParse(group);
      group->Update= (void (*) (GROUP*, int, void *, double, double))langevin_Update; 
      group->start = (void (*) (GROUP*, int, void *, double, double))langevin_Update; 
   }
   assert(group != NULL); 
   free(string); 

   LANGEVIN_PARMS *parms = (LANGEVIN_PARMS *)(group->parm); 
   object_get(obj, "vcm", &parms->vcm, WITH_UNITS, 3, "0.0 0.0 0.0","l/t",NULL);
   testForRandomKeyword(obj); 
   SYSTEM *sys = system_getSystem(NULL); 
   parms->random=system_getRandom(sys);
   //parms->g=NULL; 
   if (parms->random == NULL) missingRandomError("Langevin thermostat group");
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
