#ifndef LCG64_H
#define LCG64_H

#include <mpi.h>
#include "gid.h"
#include "random.h"

typedef struct lcg64_parm_st
{
   LONG64 state;
   unsigned int multID, prime;
} LCG64_PARM;

void lcg64_init(RANDOM *object);
double lcg64(LCG64_PARM *parm);
TWO_VECTOR lcg64_2(LCG64_PARM *parm);
MPI_Datatype lcg64_MPIType(void);
char *lcg64_parse(char *string, LCG64_PARM *parm);
void lcg64_default(unsigned n, LONG64 seed, LCG64_PARM *parm);
int lcg64_checkValue(LCG64_PARM *parm);
char *lcg64_write(LCG64_PARM *parm);
unsigned lcg64_bread (unsigned char* in , LCG64_PARM *parm);
unsigned lcg64_bwrite(unsigned char* out, LCG64_PARM *parm);
void lcg64_new_state(LONG64 seed, LCG64_PARM *parm);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
