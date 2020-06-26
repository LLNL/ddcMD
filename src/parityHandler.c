#include "parityHandler.h"

#ifdef BGL
#include <rts.h>
#include <errno.h>
static struct interrupt intrinfo;
#endif

#include <stdio.h>
#include <assert.h>

#include "mpiUtils.h"

static int _parityStatus = 0;
static int _errorsCaught = 0;
static enum ParityMode _parityMode = PARITY_ABORT;

typedef struct handler_parms_st
{
   unsigned parity_ea;
   unsigned parity_iar;
} handler_parms;

/** The first two parameters of the function are not used but are
 *  required by the signature for interrupt handers. */
void parityhandler(unsigned int interrupts, unsigned int *enabled, void *parms)
{
   ++_errorsCaught;
   // cast parms to handler variable
   handler_parms * handler_input =  (handler_parms *) parms;  
   printf("User parity handler invoked on task %d.\n"
	  "EA =%x IAR =  %x\n",
	  getRank(0),
	  handler_input->parity_ea,
	  handler_input->parity_iar ); 
   fflush(stdout);
   switch (_parityMode)
   {
     case PARITY_ABORT:
      printf("Mode = ABORT.  Calling appParityAbort() \n");
      appParityAbort();
      assert(1 == 0); // appParityAbort should not return.
      break;
     case PARITY_EXIT:
      printf("Mode = EXIT.  Calling appParityExit() \n");
      appParityExit();
      assert(1 == 0); // appParityExit should not return.
      break;
     case PARITY_RECOVER:
      _parityStatus = 1;
      break;
     default:
      printf("Unrecognized ParityMode.  This can't happen.\n");
      assert( 1 == 0);
   }
}

int parityStatus(void)
{
   return _parityStatus;
}

void parityReset(void)
{
   _parityStatus = 0;
}

int parityErrorsCaught(void)
{
   return _errorsCaught;
}

void setParityHandler(void)
{
#ifdef BGL
   static unsigned enable_parm;
   intrinfo.group = INTERRUPT_PARITY;
   intrinfo.interrupts = 0;
   intrinfo.flags = 0;
   intrinfo.handler = parityhandler;
   intrinfo.enabled = &enable_parm;
   
   int rc = rts_interrupt_control(
      INTERRUPT_SET_ENABLE, &intrinfo, sizeof(intrinfo));

   assert(rc==0);
#endif
}

void parityMode(enum ParityMode p)
{
   _parityMode = p;
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
