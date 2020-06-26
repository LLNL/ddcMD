#ifndef EAM_FS_H
#define EAM_FS_H

#include "eam.h"
#include "potential.h"
#define CMAX 72 
typedef  struct fs_pass_parms_st { double v[CMAX],d[CMAX],r2_expansion,a,b,c,m,n,l,ro,x;} FS_PASS_PARMS ;
void eam_fs_parms(POTENTIAL *object, EAM_PARMS *parms);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
