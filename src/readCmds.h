#ifndef READ_CMDS_H
#define READ_CMDS_H

#include "ddc.h"
#include "simulate.h"

/** Support for the ddcMD_CMDS facility that is used to pass commands to
 * a running simulation through a magic file.  */

enum READ_CMD_ENUMS {CHECKPOINT=1, STOP=2, DUMP_PROFILE=4, NEW_OBJECT=8, HPM_PRINT=16, DO_ANALYSIS=32};

int  readCMDS(char *filename); 

void object_rescan(DDC* ddc, SIMULATE* simulate);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
