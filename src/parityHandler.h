#ifndef PARITYHANDLER_H
#define PARITYHANDLER_H

/** Functions related to catching parity errors on BGL.  All these
 *  functions are safe to call on any platform.
 *
 * IMPORTANT *  IMPORTANT * IMPORTANT * IMPORTANT * IMPORTANT *      
 *
 * You must turn on the BGL_APP_L1_WRITEBACK_RECOVERY environment
 * variable for the parity error handler to be invoked.  If you
 * register a handler without setting this environment variable, the
 * handler will successfully register but will never be invoked.  The
 * environment variable is not enabled by default because clean cache
 * line recover does take time to perform, and thus has performance
 * implications that many users will not want to incur.
 *
 * IMPORTANT *  IMPORTANT * IMPORTANT * IMPORTANT * IMPORTANT *      
 *
 * The application is expected to supply definitions for the functions
 * void appParityAbort() and void appParityExit().  These functions
 * should not return.  See parityMode for the conditions under which
 * these functions are called.
 *
 * parityStatus:        Returns non-zero if there has been a parity error
 *                      on this task since the last call to parityReset()
 *
 * parityReset:         Resets the parity error status to zero.
 *
 * parityErrorsCaught:  Returns the number of parity errors caught on
 *                      this task.
 *
 * setParityHandler:    Sets the rts_interrupt_control to point to our
 *                      error handler (parityhandler).  Aborts 
 *                      if there is a problem setting up the interrupt.
 *                      Always returns on non-BGL hardware.
 *                 
 * parityMode:          Sets the operating mode of the parity handler.
 *                      Three modes are defined:
 *                      PARITY_ABORT:   The handler calls appParityAbort() 
 *                      PARITY_EXIT:    The handler calls appParityExit()
 *                      PARITY_RECOVER: The handler causes parityStatus()
 *                                      to return non-zero until parityReset()
 *                                      is called.
 */

enum ParityMode { PARITY_ABORT, PARITY_EXIT, PARITY_RECOVER };

/** The application must supply definitions for the following two
 *  functions! */
void appParityAbort(void);
void appParityExit(void);


int parityStatus(void);
void parityReset(void);
int parityErrorsCaught(void);
void setParityHandler(void);
void parityMode(enum ParityMode);

#endif // #ifndef PARITYHANDLER_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
