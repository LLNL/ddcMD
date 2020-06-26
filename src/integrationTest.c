#include "integrationTest.h"
#include <mpi.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "object.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "energyInfo.h"
#include "units.h"
#include "printinfo.h"
#include "state.h"
#include "testingUtils.h"

#define INTCHECK(x)  

static INTEGRATIONTEST *integrationTest = NULL;
static FILE *feng =NULL;
void zeroRpairs(SYSTEM*sys);
/*   function generated warning messages and it not used. Commented-out to avoid warning messages. 
dump_forces(SIMULATE * sim)
{
    STATE * state = sim->system->collection->state; 
    double * fx = state->fx; 
    double * fy = state->fy;
    double * fz = state->fz;  



}
*/

INTEGRATIONTEST* integrationTest_init(void *parent, char *name, MPI_Comm comm)
{
   if (name == NULL && object_exists("integrationTest", "INTEGRATIONTEST")) name = "integrationTest"; 
   if (name == NULL)  error_action("No name provided for INTEGRATION TEST object", ERROR_IN("integrationTest_init", ABORT));
   if (!object_exists(name, "INTEGRATIONTEST"))  error_action("An INTEGRATIONTEST object with name ",name,", can not be found", ERROR_IN("integrationTest_init", ABORT));
   
   //static INTEGRATIONTEST *integrationTest;
   integrationTest = (INTEGRATIONTEST *) object_initialize(name, "INTEGRATIONTEST", sizeof(INTEGRATIONTEST));
   integrationTest->parent = parent; 
   integrationTest->comm = comm;
   OBJECT* obj = (OBJECT*) integrationTest;   
   INTEGRATIONTEST *it = integrationTest;


   //whether or note positions should get printed ever iteration
   it->dump_positions = object_get(obj, "dump_positions", &(it->dump_positions)  , INT, 16, "0");

   //get list of potentials to test. should be list of pairs 
   it->nPotentialPotentialTests = object_getv(obj, "testPotentialPotential", (void *)&(it->PotentialPotentialTests) , STRING, IGNORE_IF_NOT_FOUND);
   for (int i =0 ; i<it->nPotentialPotentialTests ;i ++)
   {
      printf("potential name %s \n", it->PotentialPotentialTests[i]);
      if (!object_exists(it->PotentialPotentialTests[i], "POTENTIAL"))
      {
          potential_init(NULL,it->PotentialPotentialTests[i]);
      }
   }
   //get list of potentials to test against cached correct data. should also be list of pairs
   it->nPotentialDataTests = object_getv(obj, "testPotentialData", (void *)&it->PotentialDataTests , STRING, IGNORE_IF_NOT_FOUND);    

   //if either lists not a list of pairs, complain
   if (it->nPotentialPotentialTests%2 || it->nPotentialDataTests%2)
   {
      error_action("list of potentials outputs to test is not even", ERROR_IN("integrationTest_init", ABORT));
   }

   if (it->nPotentialPotentialTests + it->nPotentialDataTests ==0)
   {
      error_action("no tests/potentials listed", ERROR_IN("integrationTest_init", ABORT));
   }

   //simply divide nTest by two to get number of pairs
   it->nPotentialPotentialTests = it->nPotentialPotentialTests/2;
   it->nPotentialDataTests = it->nPotentialDataTests/2;

   printf("potpot %i potdat %i\n",integrationTest->nPotentialPotentialTests,integrationTest->nPotentialDataTests); 
   //printf ("integration test init has worked!!"); 
   return integrationTest;
}
/*
void reInitPotentials(SYSTEM * system)
{
   INTEGRATIONTEST *it = integrationTest;

   int nPotentials =it->nPotentialPotentialTests +it->PotentialPotentialTests;

   POTENTIAL ** iPotentials = ddcCalloc(nPotentials, sizeof(POTENTIAL *));

   for (int i =0 ;i <system->npotential; i++)
   {
      iPotentials[i] = system->potential[i]; 
   }

   int potcount =system->npotential;
   for (int i =0 ; i<2*it->nPotentialPotentialTests ;i ++) 
   {   
      printf("looking for%s \n", it->PotentialPotentialTests[i]);
      if (!object_instance_exists(it->PotentialPotentialTests[i], "POTENTIAL"))
      { 
          printf("making for%s \n", it->PotentialPotentialTests[i]);  
          POTENTIAL * p = potential_init(NULL,it->PotentialPotentialTests[i]);
          iPotentials[potcount] = potential_init(system, it->PotentialPotentialTests[i]);
          potcount++;
      }   
   }     

   for (int i =0 ; i<2*it->nPotentialDataTests ;i+=2)
   {
      if (!object_instance_exists(it->PotentialDataTests[i], "POTENTIAL"))
      {
          POTENTIAL * p = potential_init(NULL,it->PotentialDataTests[i]);
          iPotentials[potcount] = potential_init(system, it->PotentialDataTests[i]);
          potcount++;
      }
   } 

   system->potential = iPotentials;   
   char * name;
   object_get((OBJECT *) system, "neighbor", &(name), STRING, 1, NULL);
   system->neighbor = neighbor_init(system,name);
   object_get((OBJECT *) system, "volume_per_atom", &system->energyInfo.vol_per_atom, DOUBLE, 1, "-1");
}
*/

int runTests()
{     
   INTEGRATIONTEST *it = integrationTest; 
   SIMULATE * simulate1; 
   int noError=1; 
#if 0 
   int nTests =2*(it->nPotentialPotentialTests+it->nPotentialDataTests);
   int nt =1;   
   ETYPE * e1 = (ETYPE *) ddcMalloc(nTests*sizeof(ETYPE));
   ETYPE * e2 = (ETYPE *) ddcMalloc(nTests*sizeof(ETYPE));
   zeroEType(e2);
   zeroEType(e1);

   feng = fopen("feng.csv", "w");
   for (int i =0 ;i < it->nPotentialPotentialTests; i++)
   {
 		int potidx = 2*i; //index of first potential in pair
		//after the first run, we need to reset the simulation before each new run
		object_list_purge(); 
		resetStaticSpeciesVars();
		simulate1 = simulate_init(NULL,"simulate", it->comm); 

		//retrieve simulate and potential objects
		POTENTIAL * potential1 =(POTENTIAL *) object_longFind(it->PotentialPotentialTests[potidx], "POTENTIAL");
		//printf("pot name %s\n", potential1->name);

		fprintf(feng, "potential: %s\n", potential1->name);
		//TODO: this only tests one time step
		//We may need to somehow put an integrator here
		SYSTEM * sys = simulate1->system;
		simulate1->system->potential[0]->call_fsumX =1;
		zeroRpairs(sys);
		setupFirstEnergyCall(simulate1);
		setupEnergyCall(simulate1);
		char *potname1 = potential1->name;
		nt=1;
	  	for (int t = 0 ; t <nt; t++)
      {
		
		//execute first potential in pair
			printf("pre potential data\n");
			//dump_state(simulate1->system->collection->state);
         printf("fsum %i \n", potential1->call_fsumX);
			if (sys->potential[0]->call_fsumX) ddcUpdateAll(simulate1->ddc, sys, &(sys->energyInfo), 1);
			potential1->eval_potential((void *)sys, (void *)potential1->parms, &(sys->energyInfo)); 
         //if (potential1->call_fsumX) fsumX(sys, &(sys->energyInfo.virial));    
			finishEnergyCall(simulate1); 
			printf("post potential data\n");   
			//dump_state(simulate1->system->collection->state); 

			memcpy(&(e1[i]), &(sys->energyInfo), sizeof(ENERGYINFO));    
			if (getRank(0) == 0)
			{   
				printinfo(simulate1, &sys->energyInfo);
				printinfoAll(simulate1, &sys->energyInfo);
			}   

			//TODO: we need an object_remove method, so we don't keep leaking memory in object database
			//now do the same thing for second potential in pair
			object_list_purge();
			resetStaticSpeciesVars();
			SIMULATE * simulate2 = simulate_init(NULL,"simulate", it->comm);
			POTENTIAL * potential2 =(POTENTIAL *) object_longFind(it->PotentialPotentialTests[potidx+1], "POTENTIAL");
			//dump_state(simulate2->system->collection->state);
			setupFirstEnergyCall(simulate2);
			setupEnergyCall(simulate2);
			SYSTEM* sys2 = simulate2->system;
			char *potname2 = potential2->name;
			zeroRpairs(sys2);
					
			fprintf(feng, "potential: %s\n", potential2->name);
			potential2->eval_potential((void *)sys2, (void *)potential2->parms, &(sys2->energyInfo));
         //if (sys2->potential[1]->call_fsumX) fsumX(sys2, &(sys2->energyInfo.virial));
 			finishEnergyCall(simulate2);
			memcpy(&(e2[i]), &(sys2->energyInfo), sizeof(ENERGYINFO));
			
			if (getRank(0) == 0)
			{   
				printinfo(simulate2, &sys2->energyInfo);
				printinfoAll(simulate2, &sys2->energyInfo);
			} 

			UNITS *units=simulate2->printinfo->units; 
			//test correctness
			//noError =noError*compareEnergies(&e1[i],&e2[i], potname1, potname2,units, sys->nglobal,  1e-5 ) ;   
         noError =noError*compareForces(sys->collection->state,sys2->collection->state, potname1, potname2,units, sys->nglobal,  1e-3 ) ;
		} 
		//todo test PotentialData
	}
#endif 
   return noError;
   
}

double relative_error(double a, double b)
 {
   if (fabs(a)>=1e-10){
      INTCHECK(printf("a %f  b %f r%f \n", a, b,fabs((a-b)/a));)
   	return fabs((a-b)/a);
   }
   else
      INTCHECK(printf("a %f  b %f r%f \n", a, b,fabs((a-b)/a));)
		return fabs(a-b);
}

int isEqualDouble(double a, double b, double tolerance)
{
   return (relative_error(a,b)<tolerance);
}
int compareForces(STATE *s1, STATE *s2, char *potname1, char *potname2, UNITS *units,int nglobal, double tolerance )
{
   int noError = 1;
	for (int i =0 ; i<s1->nlocal; i++){
			noError = noError & isEqualDouble(s1->fx[i], s2->fx[i], tolerance);
		   fprintf(feng, "x %i : %f %f \n", i, s1->fx[i],s2->fx[i]);

			noError = noError & isEqualDouble(s1->fy[i], s2->fy[i], tolerance);
		   fprintf(feng, "y %i : %f %f \n", i, s1->fy[i],s2->fy[i]);

			noError = noError & isEqualDouble(s1->fz[i], s2->fz[i], tolerance);
		   fprintf(feng, "z %i : %f %f \n\n", i, s1->fz[i],s2->fz[i]);
	}	
   return noError;
}

//int compareEnergies(ETYPE *energyInfo1, ETYPE * energyInfo2, char *potname1, char *potname2, UNITS *units,int nglobal, double tolerance );
int compareEnergies(ETYPE *energyInfo1, ETYPE *energyInfo2, char *potname1, char *potname2, UNITS *units,int nglobal, double tolerance )
{

   THREE_SMATRIX sion;
   double Eref= 0;
   //if (getRank(0) != 0) return;
   //double eBath = energyInfo1->eBath*(units[PENERGY].convert/nglobal);
   double ekinetic = units[PENERGY].convert*(energyInfo1->rk/nglobal);
   double etot = units[PENERGY].convert*((energyInfo1->eion + energyInfo1->rk+energyInfo1->eBath)/nglobal);
   double epot = units[PENERGY].convert*(energyInfo1->eion/nglobal - Eref);
   //double temperature = units[PTEMPERATURE].convert*energyInfo1->temperature;    // temperature is unused. 
   //double pressure = units[PPRESSURE].convert*energyInfo1->pion;                // pressure is unused. 
   //voln = units[PVOLUME].convert*sys->box->volume/nglobal;
   sion = energyInfo1->sion; 
   SMATSCALE(sion, units[PPRESSURE].convert);

   double etot2, epot2, ekinetic2; //temperature2, voln2, pressure2;
   printf("a %s b %s\n", potname1, potname2);
   ekinetic2 = units[PENERGY].convert*(energyInfo2->rk/nglobal);
   etot2 =  units[PENERGY].convert*((energyInfo2->eion + energyInfo2->rk+energyInfo2->eBath)/nglobal);
   epot2 = units[PENERGY].convert*(energyInfo2->eion/nglobal - Eref);
   printf("epot %f epot2 %f \n", epot, epot2);
   printf("eion %f eion2 %f \n", energyInfo1->eion, energyInfo2->eion);
   printf("erk %f rk2 %f \n", energyInfo1->rk, energyInfo2->rk);
   printf("eBath %f eBath2 %f \n", energyInfo1->eBath, energyInfo2->eBath);
   int NoError =1; 
   if (!isEqualDouble(ekinetic, ekinetic2, tolerance))
   {
      printf("ERROR: ekinetic %s %s %f %f \n", potname1, potname2, ekinetic, ekinetic2); 
      NoError =0; 
   }
   if (!isEqualDouble(etot, etot2, tolerance))
   {   
      printf("ERROR: etotal  %s %s %f %f \n", potname1, potname2, etot, etot2); 
      NoError= 0; 
   }  
   if (!isEqualDouble(epot, epot2, tolerance))
   {   
      printf("ERROR: potential %s %s %f %f \n", potname1, potname2, epot, epot2); 
      NoError= 0;
   }  
  
   return NoError; 
   }
