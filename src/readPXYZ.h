#ifndef READ_PXYZ_H
#define READ_PXYZ_H

#include "ddc.h"

#ifdef __cplusplus
extern "C" {
#endif

// You might think it would make sense to keep the readPXYZ code next to
// the writePXYZ code.  For awhile that is what we did.  However, the
// current (8/2011) implementation of the writing code is coupled to
// SIMULATE, while the read is not.  Since the read code is called from
// ddc.c (while the write isn't) keeping the read separate allows us to
// keep ddc.c decoupled from SIMULATE.

int readPXYZ(DDC* ddc);


#ifdef __cplusplus
}
#endif

#endif
