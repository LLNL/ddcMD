#include "masters.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>		/*mallinfo call  MEMORY DEBUG */
#include "object.h"
#include "ddc.h"
#include "utilities.h"
#include "saveState.h"
#include "parityHandler.h"
#include "bglParity.h"
#include "ddcMalloc.h"
#include "io.h"
#include "transform.h"
#include "ddcenergy.h"
#include "heap.h"
#include "units.h"
#include "codata.h"
#include "ioUtils.h"
#include "readCmds.h"
#include "commandLineOptions.h"
#include "ptiming.h"
#include "external.h"
#include "hpmWrapper.h"
#include "mpiUtils.h"
#include "printinfo.h"
#include "state.h"
#include "orderSH.h"
#include "graph.h"
#include "thermalize.h"
#include "routineManager.h"
#include "three_algebra.h"
#include "integrationTest.h"
#include "unitTests.h"
#include "HAVEGPU.h"

static void thermalizeAtoms(SIMULATE* simulate, double temperature);
void adjustBox(SIMULATE *simulate);
void adjustBox1(SIMULATE *simulate);
PARTICLESET *getParticleSet(); 
void zeroAll(SYSTEM *); 
void kinetic_terms(SYSTEM *, int);
void forcetest(DDC *ddc, SYSTEM *sys);

void checkpointSimulate(SIMULATE *simulate)
{
   CreateSnapshotdir(simulate, NULL);
   simulate->system->collection->state->nlocal=simulate->system->collection->size; 
   writeRestart(simulate,0);
}

void transformMaster(void *parms_, MPI_Comm SimulateComm)
{
   DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);

   if (simulate->ntransform < 1 && getRank(0) == 0)
   {
      printf("ERROR:  No TRANSFORM objects specified in the SIMULATE object\n"); 
   }

   atStartThenExitTransforms(simulate->ntransform, simulate->transform);
   return;
}
void eightFoldMaster(void *parms_, MPI_Comm SimulateComm)
{
   DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   writeRestart8(simulate);
   return;
}
void thermalizeMaster(void *parms_, MPI_Comm SimulateComm)
{
   THERMALIZEMASTERPARMS *parms = (THERMALIZEMASTERPARMS *)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   thermalizeAtoms(simulate, parms->thermalTemp);
   return;
}
void analysisMaster(void *parms_, MPI_Comm SimulateComm)
{
   ANALYSISMASTERPARMS *parms = (ANALYSISMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   adjustBox(simulate);
   if (parms->energyCalc) firstEnergyCall(simulate); 

   CreateSnapshotdir(simulate, NULL);
   for (int i=0;i<simulate->nanalysis;i++) 
   {
      simulate->analysis[i]->eval(simulate->analysis[i]);
      simulate->analysis[i]->output(simulate->analysis[i]);
   }
   return;
}
void readWriteMaster(void *parms_, MPI_Comm SimulateComm)
{
   DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   adjustBox(simulate);
   firstEnergyCall(simulate); 
   CreateSnapshotdir(simulate, NULL);
   if (TEST0(simulate->loop,simulate->checkpointrate ))
 
  {
      writeRestart(simulate,1);
   }
   if (TEST0(simulate->loop,simulate->snapshotrate ))  
   {
      writeBXYZ(simulate);
      writeqlocal(simulate); 
      writePXYZ(simulate->ddc, simulate);
   }
   for (int i=0;i<simulate->nanalysis;i++) 
   {
      simulate->analysis[i]->eval(simulate->analysis[i]);
      simulate->analysis[i]->output(simulate->analysis[i]);
   }
   return;
}
void testForceMaster(void *parms_, MPI_Comm SimulateComm)
{
   DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   SYSTEM *sys=simulate->system; 
   adjustBox(simulate);
   firstEnergyCall(simulate); 
   forcetest(simulate->ddc,sys); 
}
void testPressureMaster(void *parms_, MPI_Comm SimulateComm)
{
   DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS*)parms_; 
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   SYSTEM *sys=simulate->system; 
   STATE *state = sys->collection->state; 
   for (int i=0;i<state->nlocal;i++) state->vx[i]=state->vy[i]=state->vz[i]=0.0; 
   double vol0 = box_get_volume(NULL);
   THREE_MATRIX h0 = box_get_h(NULL);
   adjustBox(simulate);
   firstEnergyCall(simulate); 
   double cV = units_convert(1.0, NULL,"Ang^3")/sys->nglobal;
   double cP = units_convert(1.0, NULL,"eV/Ang^3");
   double cE = units_convert(1.0, NULL,"eV")/sys->nglobal;

   sys->pCalculate=7; 
   int nn = 1;
   double pressure[2*nn+1],energy[2*nn+1],volume[2*nn+1]; 
   FILE *file[3]={NULL,NULL,NULL}; 
   if (getRank(0) == 0) 
   {
      file[0] = fopen("pressure0.data","w"); 
      file[1] = fopen("pressure1.data","w"); 
      file[2] = fopen("pressure2.data","w"); 
   }
   for (int dir=0;dir<3;dir++)
   {
      double delta = 0.4;
      for (int kk=0;kk< 12; kk++)
      {
         fprintf(file[dir],"#delta =  %f\n",delta); 
         for (int i=0;i<=2*nn+1;i++) 
         {
            volume[i]  =  vol0 + (delta*(i-nn))/cV ;
            double lambda= volume[i]/vol0; 
            THREE_MATRIX h = h0; 
            if (dir==0) h.xx = h0.xx*lambda; 
            if (dir==1) h.yy = h0.yy*lambda; 
            if (dir==2) h.zz = h0.zz*lambda; 
            box_put(NULL, HO, &h);
            THREE_MATRIX h1; 
            h1=box_get_h(NULL);

            //box_put(NULL,VOLUME,volume+i); 
            adjustBox1(simulate);
            ddcenergy(simulate->ddc, sys, 1);
            ETYPE energyInfo = sys->energyInfo; 
            //pressure[i] = -TRACE(energyInfo.sion)/3.0;   //sion OK
            if (dir==0) pressure[i] = -energyInfo.sion.xx;   //sion OK
            if (dir==1) pressure[i] = -energyInfo.sion.yy;   //sion OK
            if (dir==2) pressure[i] = -energyInfo.sion.zz;   //sion OK
            energy[i] = energyInfo.eion; 
         }
         if (getRank(0) == 0) 
         {
            for (int i=1;i<2*nn;i++) 
            {
               double pressureNumerical = -(energy[i+1]-energy[i-1])/(volume[i+1]-volume[i-1]); 
               double err = pressure[i]-pressureNumerical; 
               fprintf(file[dir],"%25.15e",volume[i]*cV); fflush(file[dir]); 
               fprintf(file[dir]," %25.15e %25.15e %25.15e %e %e",energy[i]*cE,pressure[i]*cP,pressureNumerical*cP,err*cP,err/pressure[i]); fflush(file[dir]); 
               fprintf(file[dir],"\n"); 
            }
         }
         delta *= 0.5;
      }
      fclose(file[dir]); 
   }
}

void integrationTestMaster(void *parms_, MPI_Comm SimulateComm)
{
   SIMULATEMASTERPARMS *parms = (SIMULATEMASTERPARMS *)parms_;
   SIMULATE* simulate = simulate_init(NULL,parms->common.simulateName,SimulateComm);
   adjustBox(simulate);
   //INTEGRATIONTEST * intTest = integrationTest_init(NULL,NULL,SimulateComm);
   integrationTest_init(NULL,NULL,SimulateComm);
   int success = runTests();
   if (success)
   {
      printf("All integration tests passed\n"); 
   } 
   else
   {
      printf("Integration tests failed\n");
      exit(1);
   }

   //we want to selectively call potentials based on what's listed in the integration test master

   //check if there is a TEST object in the object.data, if not, complain and exit
   //put these in integrationtest.c
   //if !(object_exists("integrationTest ","INTEGRATIONTEST")) error_action("No name provided for INTEGRATIONTEST object", ERROR_IN("integrationtestmaster", ABORT));
   //static INTEGRATIONTEST *integrationTest;
   //integrationTest = (INTEGRATIONTEST *) object_initialize(name, "INTEGRATIONTEST", sizeof(INTEGRATIONTEST));
   //the test object has a tests keyword which points to a list of test objects, get this 

   //complain if that test object doesn't exist

   //loop over every test object in the tests list
   //for (int i=0;i< nTests; i++)
   //{

   //each test object should have some sort of runTest method; if this was c++ it'd be an abstract class

   //each test should have either two  potentials or potential and data file

   //call potential one, store the energies

   //call potential two or read in energies from a cached  data file, store energies

   //print result of test
   //}

   //count good tests and bad test, say which were good and which were bad, dump a file maybe? 
}

void unitTestMaster(void * parms, MPI_Comm SimulateComm)
{
   //SIMULATEMASTERPARMS *smParms = (SIMULATEMASTERPARMS *)parms;
   //SIMULATE* simulate = simulate_init(NULL,smParms->simulateName,SimulateComm);
   //adjustBox(simulate);
   //UNITTEST * uTest = unitTest_init(NULL,NULL,SimulateComm);
   printf("running u tests\n");
   RunAllUnitTests();
}


static FILE* ldblfile=NULL;
int64_t   findEndLoop(SIMULATE *simulate)
{
   unsigned  endloop = simulate->maxloop; 
   for (int64_t loop=simulate->loop+1;loop<simulate->maxloop; loop++) 
   {
      for (int ii=0;ii<simulate->nanalysis;ii++) 
      {
         if (TEST0(loop, simulate->analysis[ii]->eval_rate)) endloop = loop;
         if (TEST0(loop, simulate->analysis[ii]->outputrate))endloop = loop;
      }
      for (int ii=0; ii<simulate->ntransform;++ii) if (TEST0(loop, simulate->transform[ii]->rate) ) endloop =loop; 
      if (TEST0(loop, simulate->printrate))  endloop = loop;
      if (TEST0(loop, simulate->snapshotrate))  endloop = loop;
      if (TEST0(loop, simulate->checkpointrate)) endloop = loop; 
      if((simulate->dataBroker != NULL) && TEST0(loop,  simulate->dataBroker->updateRate)) endloop = loop; 
      if (endloop != simulate->maxloop) break; 
   }
   return endloop; 
}
void  doLoopProlog(SIMULATE *simulate, int flag) 
{
   SYSTEM *sys=simulate->system;
   STATE *state = sys->collection->state;
   HPM_Start("LOOP");
   profile(LOOP, START);
   profileN(sys->nlocal);
   profileN2Search(sys->neighbor->nSearch);
   profileN2(sys->neighbor->npairs);
   profileN2R(sys->neighbor->nRpair);
   state->nlocal = sys->nlocal;
   state->nion = sys->nion;
}
void  doAnalysis(SIMULATE *simulate, int flag) 
{
   for (int i=0;i<simulate->nanalysis;i++) 
   {
      if (TEST0(simulate->loop, simulate->analysis[i]->eval_rate) || (flag & DO_ANALYSIS) ) simulate->analysis[i]->eval(simulate->analysis[i]);
      if (TEST0(simulate->loop, simulate->analysis[i]->outputrate)|| (flag & DO_ANALYSIS) ) simulate->analysis[i]->output(simulate->analysis[i]);
   }
}
void  doLoadBalance(SIMULATE *simulate, int flag)
{
   if ( TEST0(simulate->loop, simulate->ddc->loadBalance->rate) )
   {
      // Only open the ldbl file if we might write something to it.
      // This isn't the perfect solution as not all load balance
      // methods write, but at least this way we don't open it for
      // every simulation.
      if (getRank(0)==0) printf("load balance at loop =%"PRId64"\n",simulate->loop);
      if (ldblfile == NULL && getRank(0) == (unsigned)0) ldblfile = fopen("ldbl", "a");
      simulate->ddc->loadBalance->writeFunction( simulate->ddc->loadBalance, simulate->system->loop, simulate->time, ldblfile);
   }
}
void doCheckpoint(SIMULATE *simulate, int flag) 
{
   if (TEST0(simulate->loop, simulate->checkpointrate) || flag & CHECKPOINT ) 
   {
      CreateSnapshotdir(simulate, NULL);
      writeRestart(simulate,1);
      mem_debug("writeRestart");
   }
}
void doDumpProfile(SIMULATE *simulate,int flag)       
{
   if (TEST0(simulate->loop, simulate->snapshotrate) || TEST0(simulate->loop, simulate->checkpointrate) || (flag & DUMP_PROFILE))
   {
      CreateSnapshotdir(simulate, NULL);
      int dataLength = sizeof(double) + strlen(simulate->snapshotdir) +8+1; 
      char data[dataLength]; 
      char *profileFilename = data+sizeof(double); 
      double eloop=profileGetLoopTime(); 
      copyBytes(data, &eloop, sizeof(double)); 
      sprintf(profileFilename, "%s/profile",simulate->snapshotdir);
      callRoutine(NULL,slave_dumpprofile,dataLength,data); 
      dumpprofile(profileFilename);
   }
}
void doSnapshot(SIMULATE *simulate, int flag)
{
   if ( simulate->startLoop < simulate->loop)
   {
      if (TEST0(simulate->loop, (SIGNED64)simulate->snapshotrate)) 
      {
         CreateSnapshotdir(simulate, NULL);
         writeBXYZ(simulate);
         writeqlocal(simulate); 
         writePXYZ(simulate->ddc, simulate);
      }  
   }
}
void    doDataBroker(SIMULATE *simulate,int flag)
{
   if(simulate->dataBroker != NULL  && simulate->dataBroker->useDb)
   {
      if(TEST0(simulate->loop, simulate->dataBroker->updateRate))
      {
         int typeflag;
         typeflag = strcmp(simulate->dataBroker->type, "DB") == 0;
         typeflag |= strcmp(simulate->dataBroker->type, "TRANS_CHARMM") == 0;  
         if(typeflag) postTotalEnergies(simulate->dataBroker, simulate->system->energyInfo, simulate->loop);
         typeflag = strcmp(simulate->dataBroker->type, "HYCOP") == 0;
         if(!typeflag) updateDataBrokerContent(simulate->dataBroker, simulate->system, simulate->loop, 1);
      }
   }
}

void simulateMaster(void *parms, MPI_Comm SimulateComm)
{
   SIMULATEMASTERPARMS *smParms = (SIMULATEMASTERPARMS *)parms;
   char *simulateName = smParms->common.simulateName; 
   SIMULATE* simulate = simulate_init(NULL,simulateName,SimulateComm);
   SYSTEM *sys = simulate->system; 
   //STATE *state = sys->collection->state; 
   adjustBox(simulate);
   if(simulate->dataBroker != NULL)
   {
      int flag = strcmp(simulate->dataBroker->type, "HYCOP") == 0;
      if(simulate->dataBroker->useDb && !flag)
      {
         int cdberr = connectToDB(simulate->dataBroker);
         printf("Connecting to Data Broker\n");
         if(cdberr>0)
            printf("ERROR CONNECTING TO DB: %d\n", cdberr);
      }
   }

   int gpu_integrate = simulate->integrator->uses_gpu;
   if(gpu_integrate)
   {
      GPUCODE(sendGPUState(sys, sys->nion);)
      GPUCODE(sendForceVelocityToGPU(sys, sys->nlocal);)
   }

   firstEnergyCall(simulate); 
   if(gpu_integrate==1)
   {
      GPUCODE(sendForceEnergyToHost(sys,&sys->energyInfo);)
      GPUCODE(sendForceVelocityToHost(sys, sys->nlocal);)
      kinetic_terms(sys,1);
      eval_energyInfo(sys); 
   }
   //   if (getRank(0) ==0) printf("EION = %f \n",sys->energyInfo.eion); 
   profile(NOFIELD, ZERO);  //zero timers. 
   FILE* graphfile=NULL;
   if (getRank(0) == 0 && simulate->printinfo->printGraphs != 0) graphfile = fopen("graphs", "a");
   graphWrite(simulate, graphfile);
   FILE* ldblfile=NULL; 


   if (getRank(0) == 0) 
   {
      printinfo(simulate, &sys->energyInfo);
      printinfoAll(simulate, &sys->energyInfo);
   }

   timestamp("Starting MD loops");
   time_t current_time= time(NULL); 
   //timestamp("calling saveState()...\n");
   //saveState(simulate);
   parityMode(PARITY_RECOVER);
   for (int i=0;i<simulate->nanalysis;i++) simulate->analysis[i]->startup(simulate->analysis[i]);

   simulate->startLoop = simulate->loop; 
   for (; simulate->loop < simulate->maxloop;) /*  main MD loop  */
   {		
      doLoopProlog(simulate,0); 
      int64_t endLoop = findEndLoop(simulate);

      if(simulate->dataBroker != NULL)
      { 
         int flag = strcmp(simulate->dataBroker->type, "HYCOP") == 0;
         if(simulate->dataBroker->useDb && TEST0(simulate->loop, simulate->dataBroker->updateRate) && !flag)
         {
            updateDataBrokerContent(simulate->dataBroker, simulate->system, simulate->loop, 1);
         }
      }

      HPM_Start("Integrator");

      for (; simulate->loop < endLoop;) 
      {
         profile(MDSTEP, START);                          //  All Variables at t ;  
         if (simulate->itype == MD) simulate->integrator->eval_integrator(simulate->ddc, simulate, simulate->integrator->parms);
         profile(MDSTEP, END);                          //  All Variables at t = t+h ;  
      }
      if(gpu_integrate==1)
      {
         GPUCODE(sendForceEnergyToHost(sys,&sys->energyInfo);)
         GPUCODE(sendForceVelocityToHost(sys, sys->nlocal);)
         GPUCODE(sendPosnToHost(sys, sys->nlocal); )
      }
      kinetic_terms(sys, 1);
      eval_energyInfo(sys);

      atRateTransforms(simulate->ntransform, simulate->transform);

      HPM_Stop("Integrator");
      HPM_Stop("LOOP");
      callRoutine(NULL, slave_parityFailure, 0, NULL);
      if (parityFailure(simulate, simulate->ddc) != 0)
      {
         profile(LOOP, END);
         continue;
      }
      HPM_Start("loopio");
      profile(LOOPIO, START);
      int flag = 0;
      if ( !isfinite(sys->energyInfo.eion)) 
      {
         flag= STOP; 
         printf("eion = %e is bad. Simulation is being killed at loop = %"PRId64"\n",sys->energyInfo.eion,simulate->loop);
         printinfoAll(simulate, &sys->energyInfo);
      }
      if (flag != STOP) 
      {
         if (TEST0(simulate->loop, simulate->printrate))
         {

            mem_debug("potential");
            printinfoAll(simulate, &sys->energyInfo);
            //if (getRank(0) == 0) ddcMemSummary(stdout);
            graphWrite(simulate, graphfile);
            flag = readCMDS("./ddcMD_CMDS");
         }
         if (!(flag & (CHECKPOINT | STOP) ) )
         {
            if (getRank(0) == 0)
            {
               int deltaTime=120; 
               time_t last_time = current_time; 
               current_time = time(NULL); 
               deltaTime = MAX(2*(current_time-last_time), deltaTime); 
               if (smParms->stopTime > -1 && smParms->stopTime < current_time+deltaTime) 
               {
                  flag |= (CHECKPOINT | STOP); 
                  timestamp("Program is stopping to avoid exceeding time limit"); 
               }
            }
            MPI_Bcast(&flag,1,MPI_INT,0, SimulateComm); 
         }
         doCheckpoint(simulate, flag); 
         doAnalysis(simulate, flag);
         doSnapshot(simulate, flag);
         doLoadBalance(simulate, flag);
      }

      profile(LOOPIO, END);
      profile(LOOP, END);

      doDumpProfile(simulate,flag);       
      HPM_Stop("loopio");
      if (flag & HPM_PRINT)
      {
         HPM_Stop("Total");
         HPM_Print("hpm.data");
         HPM_Start("Total");
      }
      if (flag & NEW_OBJECT)
      {
         object_Bcast(0, SimulateComm);
         object_rescan(simulate->ddc, simulate);
      }
      if (flag & STOP ) break; 
      //HPM_Group_Roll(); 
      doDataBroker(simulate,flag);
   }

   timestamp("Finished MD loops");
   if (!TEST0(simulate->loop, simulate->printrate)) 
   {
      if(gpu_integrate==1)
      {
         GPUCODE(sendForceEnergyToHost(sys,&sys->energyInfo);)
         GPUCODE(sendForceVelocityToHost(sys, sys->nlocal);)
      }
      kinetic_terms(sys,1); 

      printinfo(simulate, &sys->energyInfo); 
      timestamp("Finished printinfo");
      printinfoAll(simulate, &sys->energyInfo); 
      timestamp("Finished printinfoAll");
      graphWrite(simulate, graphfile);
   }
   if ( simulate->startLoop < simulate->loop) doDumpProfile(simulate,DUMP_PROFILE);       
   if(simulate->dataBroker != NULL)
      if(simulate->dataBroker->useDb && simulate->dataBroker->stopDb) shutdownDB(simulate->dataBroker); //Close DB if we are in charge of managing the DB

   timestamp("Begin closing files");
   timestamp("  printinfo file");
   printinfo_close(simulate); 
   timestamp("  graph file");
   if (getRank(0) == 0 && simulate->printinfo->printGraphs !=0) fclose(graphfile);
   timestamp("  ldbl file");
   if (getRank(0) == 0 && ldblfile != NULL) fclose(ldblfile);
   timestamp("End closing files");
   heap_deallocate(); 
}

void thermalizeAtoms(SIMULATE* simulate, double temperature)
{
   timestamp("start retherm\n");
   simulate->loop=0; 
   simulate->time=0; 

   THERMALIZE_PARMS* parms = thermalize_init();
   parms->method = THERMALIZE_BOLTZMANN;
   parms->seed = generateRandomSeed();
   thermalize_setTemperature(parms, THERMALIZE_GLOBAL, &temperature, NULL, 1);
   thermalize(simulate->system, parms);
   thermalize_destroy(parms);
   ddcFree(parms);

   CreateSnapshotdir(simulate, "snapshotThermal");
   writeRestart(simulate,1);
   timestamp("Exit Thermalize");
}
void firstEnergyCall(SIMULATE *simulate)
{
   SYSTEM *sys = simulate->system; 
   COLLECTION *collection = sys->collection;
   STATE *state = collection->state; 
   ddc_put(simulate->ddc, PAIRMODE, LABEL);
   ddc_put(simulate->ddc, CORNER, &sys->box->corner);
   ddc_put(simulate->ddc, DDCNLOCAL, sys->nlocal);
   simulate->ddc->update = 3;

   
  if ((simulate->accelerator != NULL) && simulate->accelerator->itype == GPU_CUDA)
   {
   GPUNLIST *gnlist= sys->collection->gnlist; 
   GPUCODE(allocPages(gnlist, sys->nion, 1,simulate->ddc->rcut );)
   }
   //serializeSpecies(sys, sys->nlocal);

   resize(collection->size, 2, collection->state); // ensure that any quantities registered with particleRegisterinfo are allocated before we try to use them in ddcenergy.
   for (int kk=0; kk<sys->ngroup; kk++) 
   {
      if (sys->group[kk]->start != NULL) sys->group[kk]->start(sys->group[kk],START_GROUP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
      sys->group[kk]->Update1(sys->group[kk],START_GROUP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
   }

   if(simulate->dataBroker != NULL){
      if(simulate->dataBroker->useDb && TEST0(0, simulate->dataBroker->updateRate))
      {
         updateDataBrokerContent(simulate->dataBroker, simulate->system, simulate->loop, 0);
      }
   }

   timestamp("Starting First Call to ddcenergy");
   ddcenergy(simulate->ddc, sys, 1);
   timestamp("Finished First Call to ddcenergy");

   resize(simulate->ddc->number_particles, 2, state);
   for (int kk=0; kk<sys->ngroup; kk++) 
   {
      sys->group[kk]->Update(sys->group[kk],FRONT_TIMESTEP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
   }
}
void adjustBox(SIMULATE *simulate) 
{
   THREE_MATRIX hfac; 
   SYSTEM *sys = simulate->system; 
   box_put(sys->box,BOX_TIME,(void *)&simulate->time);
   box_get(sys->box,HFAC,(void *)&hfac);
   double *rx = sys->collection->state->rx;
   double *ry = sys->collection->state->ry;
   double *rz = sys->collection->state->rz;
   double *vx = sys->collection->state->rx;
   double *vy = sys->collection->state->ry;
   double *vz = sys->collection->state->rz;
   for (unsigned kk = 0; kk < sys->nlocal; kk++)
   {
      THREE_VECTOR u,v; 
      u.x=rx[kk]; u.y=ry[kk]; u.z=rz[kk];
      v = matrix_vector(hfac,u); 
      rx[kk]=v.x; ry[kk]=v.y; rz[kk]=v.z;
      u.x=vx[kk]; u.y=vy[kk]; u.z=vz[kk];
      v = matrix_vector(hfac,u); 
      vx[kk]=v.x; vy[kk]=v.y; vz[kk]=v.z;
   }
}
void adjustBox1(SIMULATE *simulate) 
{
   THREE_MATRIX hfac; 
   SYSTEM *sys = simulate->system; 
   box_get(sys->box,HFAC,(void *)&hfac);
   double *rx = sys->collection->state->rx;
   double *ry = sys->collection->state->ry;
   double *rz = sys->collection->state->rz;
   double *vx = sys->collection->state->vx;
   double *vy = sys->collection->state->vy;
   double *vz = sys->collection->state->vz;
   for (unsigned kk = 0; kk < sys->nlocal; kk++)
   {
      {
      THREE_VECTOR u = {rx[kk],ry[kk],rz[kk]}; 
      THREE_VECTOR v = matrix_vector(hfac,u); 
      rx[kk]=v.x; ry[kk]=v.y; rz[kk]=v.z;
      }
      {
      THREE_VECTOR u = {vx[kk],vy[kk],vz[kk]}; 
      THREE_VECTOR v = matrix_vector(hfac,u); 
      vx[kk]=v.x; vy[kk]=v.y; vz[kk]=v.z;
      }
   }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
