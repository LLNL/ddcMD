// accumulate a simple time average of the force on certain particles
// append the values to an output file
//
#include <stdlib.h>
#include <string.h>
#include <assert.h>

//#include "box.h"
#include "object.h"
#include "pio.h"
#include "io.h"
#include "crc32.h"
#include "ddcMalloc.h"
#include "units.h"
#include "mpiAlgorithm.h"
#include "pinfo.h"
#include "ioUtils.h"
#include "system.h"
#include "three_algebra.h"

#include "forceAverage.h"
#include "format.h"

 int getRank(int);
 void forceAverage_clear(ANALYSIS* analysis);

////////////////////////////////////////////////////////////////////////
FORCE_AVERAGE_PARMS *forceAverage_parms(ANALYSIS *analysis)
{
   FORCE_AVERAGE_PARMS *parms = ddcMalloc(sizeof(FORCE_AVERAGE_PARMS)); 

// SIMULATE* simulate =(SIMULATE *)analysis->parent; 
// SYSTEM* sys=simulate->system;

   OBJECT* obj = (OBJECT*) analysis;

   object_get(obj, "eval_rate", &parms->_eval_rate, INT, 1, "1");
   object_get(obj, "outputrate", &parms->_outputrate, INT, 1, "1");
   object_get(obj, "filename",   &parms->_filename,  STRING, 1, "force.dat");
   object_get(obj, "data_type",   &parms->_data_type,  STRING, 1, "force");

   //////////////////////////////////////////////////
   // Read the GID of all particles to be tracked. //
   //////////////////////////////////////////////////
   parms->_nParticles = object_getv(obj, "idParticle",
                                      (void*)&parms->_idParticle,
                                      U64, ABORT_IF_NOT_FOUND);

   /////////////////////////////////////////////////////
   // Create storage for all histograms on all tasks. //
   /////////////////////////////////////////////////////
   parms->_forbar  = ddcMalloc(parms->_nParticles*sizeof(THREE_VECTOR));

   for (int ii=0;ii<parms->_nParticles; ii++)
        VSET(parms->_forbar[ii],0.,0.,0.);

   parms->_nSnapSum = 0;

   // print header
   if (getRank(0)==0)
      {
       parms->_file = fopen(parms->_filename,"a");
  
       fprintf(parms->_file,"# loop  time  ");
       for (int i=0;i<parms->_nParticles;i++)
           fprintf(parms->_file,"%llu  ",(unsigned long long int)parms->_idParticle[i]);
       fprintf(parms->_file,"\n");
       fflush(parms->_file);
   
      }

   return parms;
}

////////////////////////////////////////////////////////////////////////
void forceAverage_eval(ANALYSIS* analysis)
{

   FORCE_AVERAGE_PARMS* parms = (FORCE_AVERAGE_PARMS *)analysis->parms; 

   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state; 
   int nlocal = simulate->system->nlocal;

   int taskNum = getRank(0);

   THREE_VECTOR sndbuf[1], recbuf[1];
   MPI_Status * MPIStatus;
   MPIStatus = (MPI_Status*) ddcMalloc(sizeof(MPI_Status));
   MPI_Request* MPIRequest;
   MPIRequest = (MPI_Request*) ddcMalloc(sizeof(MPI_Request));

   ///////////////////////////////////////
   // All tasks search nlocal particles //
   // for the iC'th particle in list.   //
   ///////////////////////////////////////
   for (int iC=0; iC<parms->_nParticles;iC++)
       {
        for (int ii=0; ii<nlocal; ii++)
            {
             ////////////////////////////////////////////////
             // One task finds the iC'th central particle! //
             ////////////////////////////////////////////////
             if (parms->_idParticle[iC] == state->label[ii])
                {
                 ////////////////////////////////////////////////
                 // Finding task: pack projectile information. //
                 // Send projectile information to task #0.    //
                 ////////////////////////////////////////////////
                 if (strcasecmp(parms->_data_type, "force") == 0)
                     {VSET(sndbuf[0],state->fx[ii],state->fy[ii],state->fz[ii]);}
                 else if (strcasecmp(parms->_data_type, "velocity") == 0)
                     {VSET(sndbuf[0],state->vx[ii],state->vy[ii],state->vz[ii]);}
                 else if (strcasecmp(parms->_data_type, "position") == 0)
                     {VSET(sndbuf[0],state->rx[ii],state->ry[ii],state->rz[ii]);}
                 else
                     {if (getRank(0)==0) printf("bad forceAverage data_type\n");exit(19);}
// mod 3/5/13 - try setting tag to unique identifier iC
                 MPI_Isend((void*)sndbuf,3,MPI_DOUBLE,0,iC,COMM_LOCAL,MPIRequest);
                }
             }

        ////////////////////////////////////////////////////////
        // Task 0: blocking receive info on iC'th projectile, //
        //         accumulate force over time.                //
        ////////////////////////////////////////////////////////
        if (taskNum==0)
           {
            MPI_Recv((void*)recbuf,3,MPI_DOUBLE,
                     MPI_ANY_SOURCE,iC,COMM_LOCAL,
                     MPIStatus);
           // parms->_forbar[iC] = vector_sadd(1.,parms->_forbar[iC],recbuf[0]);
            parms->_forbar[iC].x += recbuf[0].x;
            parms->_forbar[iC].y += recbuf[0].y;
            parms->_forbar[iC].z += recbuf[0].z;
           }

       }

 parms->_nSnapSum++;

 ddcFree(MPIStatus);
 ddcFree(MPIRequest);
 return ;
}

////////////////////////////////////////////////////////////////////////
void forceAverage_output(ANALYSIS* analysis)
{
 //////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////
 SIMULATE* simulate =(SIMULATE *)analysis->parent; 

 FORCE_AVERAGE_PARMS* parms = (FORCE_AVERAGE_PARMS *)analysis->parms; 

 if (parms->_nSnapSum*parms->_eval_rate < parms->_outputrate) return;

 //double convert;
 //if (strcasecmp(parms->_data_type, "force") == 0)
 //   {convert = units_convert(1.0,NULL,"l/t^2");}
 //else if (strcasecmp(parms->_data_type, "velocity") == 0)
 //   {convert = units_convert(1.0,NULL,"velocity");}
 //else if (strcasecmp(parms->_data_type, "position") == 0)
 //   {convert = units_convert(1.0,NULL,"l");}

 if (getRank(0)==0)
    {
     fprintf(parms->_file,loopFormat(),simulate->loop);
     fprintf(parms->_file," %10.6f ",simulate->time);
     for (int iC=0; iC<parms->_nParticles;iC++)
         {
          fprintf(parms->_file," %12.7f %12.7f %12.7f ",
//                (unsigned long long int)parms->_idParticle[iC],
                  parms->_forbar[iC].x/(double)parms->_nSnapSum,
                  parms->_forbar[iC].y/(double)parms->_nSnapSum,
                  parms->_forbar[iC].z/(double)parms->_nSnapSum);

         }
      fprintf(parms->_file,"\n");
      fflush(parms->_file);
  
     }

 ////////////////////////////////////////////////
 // Re-initialize the histograms on all tasks. //
 ////////////////////////////////////////////////
 forceAverage_clear(analysis);

}

////////////////////////////////////////////////////////////////////////
void forceAverage_clear(ANALYSIS* analysis)
{
 FORCE_AVERAGE_PARMS* parms = (FORCE_AVERAGE_PARMS *)analysis->parms;

 /////////////////////////////////
 // Re-initialize on all tasks. //
 /////////////////////////////////
 for (int ii=0;ii<parms->_nParticles; ii++)
        VSET(parms->_forbar[ii],0.,0.,0.);

 parms->_nSnapSum = 0;

}

