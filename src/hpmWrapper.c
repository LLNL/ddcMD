#include "hpmWrapper.h"

#ifdef BGL
#define HAS_HPM
// function definitions are supplied by hpm/hpm.a
#endif


// libmpihpm with block capability hangs for >20 min at MPI_Finalize
// before producing output for all blocks.  So, for now this code can't
// be used.
#ifdef BGP
#define HAS_HPM
// function definitions are supplied by hpm_bgp/libhpm.a
#endif


// Supply empty functions for platforms that don't have an hpm library.
#ifndef HAS_HPM
void HPM_Init(void) {}
void HPM_Start(char* s) {}
void HPM_Stop(char* s) {}
void HPM_Print(char* filename) {}
// getHPM_FlopCnt is not currently used in the code so I'm not going to
// provide a defintion.  Note that the definition of this function for
// BG/L in the hpm subdirectory includes an MPI_Reduce on
// MPI_COMM_WORLD.  This could be problematic depending on the context
// from which it is called.
//int getHPM_FlopCnt(char* s) {return 0;}
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
