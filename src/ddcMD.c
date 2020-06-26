//
// ddcMD is an highly scalable particle integrator
//
// ************************************************************/
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>      /*mallinfo call  MEMORY DEBUG */
#include <sys/time.h>
#include <string.h>
//#include <fenv.h>
#include "object.h"
#include "utilities.h"
#include "ioUtils.h"
#include "commandLineOptions.h"
#include "ptiming.h"
#include "routineManager.h"
#include "units.h"
#include "codata.h"
#include "masters.h"
#include "rectimer.h"
#include "primes.h"
#include "mpiUtils.h"
#include "hpmWrapper.h"
#include "parityHandler.h"
#include "hardwareInfo.h"
#include "ddcMalloc.h"
#include "expandbuffer.h"
#include "io.h"

void objectSetup(void *parms, MPI_Comm comm);
void mpiStartUp(int argc, char *argv[]);
void   parityErrorReport();
void setupForGmonFiles();

void resizePrint(PFILE*file);
void version_init(int argc,char *argv[]);
void simulate(COMMAND_LINE_OPTIONS  opt, MPI_Comm SimulateComm);
void commons_init(void);

MPI_Comm COMM_LOCAL;

/**
 *  Some notes about units:
 *
 *  It is impossible to define a unit system where all of the units are "convenient" or "customary" for atomic scale simulation.  For example, you cannot have length in Angstroms, mass in AMU, time in fs, charge in electron charge, 
 *
 *  Internal Units:  The internal units are defined such that:
 *     - length is in bohr
 *     - energy is Rydberg
 *     - time is femtosecond
 *     - charge is electron charge,
 *
 *  External Units:  The external units are chosen such that
 *     - length is Angstroms
 *     - temperature is Kelvin
 *     - time is femtosecond
 *     - charge is electron charge
 *     - mass is AMU
 * internal units: bohr, ???, fs, e-charge/fs, ???, mol, candella 
 * external units Ang, AMU, fs, e- charge/fs, Kelvin, mol, candella
*/
#if TESTEXE == 0
int main(int argc, char *argv[])
{
 //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
   mpiStartUp(argc,argv); 
   COMMAND_LINE_OPTIONS  opt = parseCommandLine(argc, argv);
   checkLimits();
   commons_init();
   prime_init(30000, getRank(-1), getSize(-1));
   units_internal(a0_MKS,Rinfhc_MKS*1e-30/(a0_MKS*a0_MKS),1e-15,e_MKS/1e-15,Rinfhc_eV/kB_eV,1.0,1.0);
   units_external(1e-10,u_MKS,1e-15,e_MKS/1e-15,1.0,1.0,1.0);
   checkpointUnits("Ang","amu","fs","e/fs","K"," ","cd");
// ddcMemSetVerbose(getRank(-1)); 
   setParityHandler();
   version_init(argc,argv);

   profile(NOFIELD, INIT);
   HPM_Init();

   profile(TOTAL, START);

   HPM_Start("Total");
   WAIT(-1);
   MASTER master=masterFactory(opt);
   objectSetup(master.parms, MPI_COMM_WORLD);
   ROUTINE *routine  = routineManager_init(NULL,"routineManager",master.fcn,master.parms);
   routine->fcn(routine->parms,routine->comm);
   routineDestroy(routine); 
   HPM_Stop("Total");
   profile(TOTAL, END);

   profile(NOFIELD, AVE);
   parityErrorReport();
   //HPM_Print("hpm.data"); //ewd comment this out to avoid data overflow
   timestamp("Waiting for all tasks to complete");
   timestampBarrier("Calling MPI_Finalize",MPI_COMM_WORLD);
   setupForGmonFiles();
   MPI_Finalize();
}   /* end of main routine */
#endif
void printNetworkInfo()
{
   const thandle t_torus = rectimer_start(NULL,"main:torus");
   if (hi_hasTorus() == 1 && getRank(-1) == 0)
   {
      int nDim = hi_nTorusDim();
      int size[nDim];
      hi_torusSize(size);
      printf("Torus size %d", size[0]);
      for (int ii=1; ii<nDim; ++ii)
         printf(" x %d", size[ii]);
      printf("\n");
   }
   rectimer_stop(t_torus);
}
void mpiStartUp(int argc, char *argv[])
{   
  /* Figure out time when we first hit main, no communication yet */

  time_t t_main_enter = time(NULL);
  struct timeval tv_main;
  gettimeofday(&tv_main,NULL);

   

   const thandle t_mpiinit = rectimer_start(NULL,"main:MPI_Init");
   MPI_Init(&argc, &argv);
   rectimer_stop(t_mpiinit);


   COMM_LOCAL = MPI_COMM_WORLD; 

   /* Print timestamp and process variance of entering to main() */
   {
     const thandle t_reduce = rectimer_start(NULL,"main:reduce time");
     double thit = tv_main.tv_sec + 1e-6*tv_main.tv_usec,tmin,tmax,tavg;
     int pid,np;
     char tstr[80];

     {
       int len;
       const char *str_ctime = ctime(&t_main_enter);
       strcpy(tstr,str_ctime);
       len = strlen(tstr);
       tstr[len-1] = '\0';
     }

     MPI_Comm_rank(COMM_LOCAL,&pid);
     MPI_Comm_size(COMM_LOCAL,&np);
     MPI_Allreduce(&thit,&tmin,1,MPI_DOUBLE,MPI_MIN,COMM_LOCAL);
     MPI_Allreduce(&thit,&tmax,1,MPI_DOUBLE,MPI_MAX,COMM_LOCAL);
     MPI_Allreduce(&thit,&tavg,1,MPI_DOUBLE,MPI_SUM,COMM_LOCAL);
     tavg /= np;
     if(pid == 0) {
       printf("%s: Process %d of %d hit main, time skew: min = %+.3fs / avg = %+.3fs / max = %+.3fs\n",
              tstr,pid,np,tmin-thit,tavg-thit,tmax-thit);
       fflush(stdout);
     }
     MPI_Barrier(COMM_LOCAL);
     rectimer_stop(t_reduce);
   }
//******************************************
#if 0 
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    if (getRank(0) == 0) 
   {
      printf("%d: PID %d on %s ready for attach\n", getRank(0),getpid(), hostname); fflush(stdout);
      sleep(30);
   }
#endif
   printNetworkInfo();
}
void   parityErrorReport()
{ //scope
   int localErr = parityErrorsCaught();
   int globalErr;
   MPI_Allreduce(&localErr, &globalErr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (getRank(-1) == 0) printf("Parity errors handled: %d\n", globalErr);
} //scope

void setupForGmonFiles()
{
#ifndef PROFILE
   return; 
#endif
   char cwd[256], path[256];
   struct stat statbuf;
   int rc, mode = 0775;
   getcwd(cwd, 255);
   sprintf(path, "%s/PG", cwd);
   if (getRank(-1) == 0)
   {
   rc = stat(path, &statbuf);
   if (rc == -1 && errno == ENOENT) rc = mkdir(path, mode);
   }
#if defined(BGL) || defined(BGP) || defined(BGQ)
   rc = chdir(path);
   return;
#endif
   WAIT(-1); 
   sprintf(path, "%s/PG/%6.6d", cwd, getRank(-1));
   rc = stat(path, &statbuf);
   if (rc == -1 && errno == ENOENT) rc = mkdir(path, mode);
   rc = chdir(path);
}

void dumpprofile(char *filename)
{
   PFILE *file;
   profile(TOTAL, END);
   profile(NOFIELD, AVE);
   file = Popen(filename, "w", MPI_COMM_WORLD);
   //WAIT(-1);
   Pprintf(file,"\n\n\n************************************************************************\n\n");
   profileWrite(file, getRank(-1));
   ddcMemReport_pio(file);
   ExpandBuffersPrint_pio(file);
   resizePrint(file);
   Pclose(file);
   profile(TOTAL, ACCUM);
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
