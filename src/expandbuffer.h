#ifndef EXPANDBUFFER_H
#define EXPANDBUFFER_H

#include "expandbuffer.h"
#include <stdio.h>
#include "pio.h"

#ifdef __cplusplus
extern "C" {
#endif

void *ExpandBuffers(void *buffer, int size, int n, int increment, char *string, ...);
void *ExpandBuffersPrint(FILE*file);
void *ExpandBuffersPrint_pio(PFILE*file);
void ExpandBuffersResize(void);
void ExpandBuffersMode(char *mode_string);
int  ExpandBuffersLastLength(void);
void  ExpandBuffersFree(void *);

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
