#include "bglParity.h"

#ifdef BGL
#include <rts.h>
#endif
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "parityHandler.h"
#include "ddc.h"
#include "ddcenergy.h"
#include "saveState.h"
#include "eam_onepass.h"
#include "mpiUtils.h"

/** This function is static so that it can't be called by unsuspecting
 * clients that fail to realize there is an MPI_Allreduce over
 * MPI_COMM_WORLD hidden inside.  When the application wishes to check
 * for a parity error it should call pairtFailure() instead (along with
 * slave_parityFailure on all slave tasks). */
static int parityError(void);

/** Provides definitions for application supplied function called
 *  by the parity handler. */
void appParityAbort(void)
{
   MPI_Abort(COMM_LOCAL, 46);
}


/** Provides definitions for application supplied function called
 *  by the parity handler. */
void appParityExit(void)
{
   writeState();
   MPI_Finalize();
   exit(0);
}


/** Checks whether a parity error has been detected.  If so, the system
 *  is restored to the last saved state and a non-zero value is returned.
 *  Otherwise the state of the system is saved and zero is returned.
 *
 *  Note that this routine is called by slave_parityFailure so it may be
 *  executed by slave tasks that have no meaningful SIMULATE or DDC.
 *
 *  LIMITATIONS: Since the saveState and restoreState functions operate
 *  on simulate (which is probably NULL on slave tasks) there isn't much
 *  to work with for the slave tasks to save and restore state.  If we
 *  develop a real need to save and restore state on the slave tasks
 *  this will need to be dealt with somehow.
 */
int parityFailure(SIMULATE* simulate, DDC* ddc)
{
#ifndef BGL
   return parityError();
#else
   
   
   int rc = 0;
   if (parityError())
   {
      do
      {
         restoreState(simulate);
         eam_onepass_reset();
      }
      while (parityError() != 0);
      if (getRank(0) == 0)
         printf("System restored to step %d to recover from parity error.\n",
                simulate->loop);
      if (ddc && simulate)
      {
         ddc->lastUpdate = NO_DDC_UPDATE;
         ddc_put(ddc, DDCNLOCAL, simulate->system->nlocal);
         ddcenergy(ddc, simulate->system, 1);
      }
      rc = 1;
   }
   else
      do
         saveState(simulate);
      while (parityError() != 0);
   
   return rc;
   
#endif
}

/** Slave tasks don't have a fully implemented or even fully thought out
 *  parity recovery strategy.  What is here is mostly pretty sketchy to
 *  serve as an outline in case we really need to do something.  I
 *  expect that we will not run many large simulations with multiple
 *  routines on BG/L platforms so I'm not going to put a great deal of
 *  time into developing a full blown parity recovery framework.
 *
 *  As far as I can recall the slave tasks don't have meaningful
 *  SIMULATE or DDC object to pass to the parityFailure function.  As of
 *  this writing they also have no saved state to restore.  So, that
 *  means that at present the slave tasks just get to raise their hand
 *  and say the result was bad so the particle tasks need to rewind.
 *
 *  Also note that this function does not return a value (though it
 *  potentially could through the argument list) so at this point there
 *  is no way for the slave tasks to know that a parity error occured.
 */
void slave_parityFailure(char* data)
{
   parityFailure(NULL, NULL);
}



/** Checks whether any task has registered a parity error.  Returns
 *  non-zero if an error is detected.
 *
 *  This routine is designed to check for parity errors on *all* tasks.
 *  That means the slave tasks in a non-homogeneous decomposition are
 *  required to check for parity errors at the same time as the master
 *  tasks.
 */
int parityError(void)
{
#ifndef BGL
   return 0;
#else
   
   profile(PARITYCHECK, START);
   rts_dcache_evict_transient();
   rts_dcache_evict_normal();
   int globalErr = 0;
   int localErr = parityStatus();
   parityReset();
   MPI_Allreduce(&localErr, &globalErr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   profile(PARITYCHECK, END);
   return globalErr;
   
#endif	
}

/* Local Variables: */
/* tab-width: 3 */
/* indent-tabs-mode: nil */
/* End: */
