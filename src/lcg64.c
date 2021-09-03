#include "lcg64.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <wordexp.h>
#include "three_algebra.h"
#include "error.h"
#include "object.h"
#include "ioUtils.h"
#include "ddcMalloc.h"
#include "primes.h"
#define TWO_M64 5.4210108624275222e-20	/* 2^(-64) */
const LONG64 INIT_SEED=0x2bc6ffff8cfe166dllu;

int getRank(int);


void lcg64_init(RANDOM *random)
{
	random->itype=LCG64;
	random->sizeOfParms  = sizeof(LCG64_PARM);
	random->parse  = (char *(*)(char *, void *))lcg64_parse;
	random->defaultValue  = (void (*)(unsigned , LONG64 , void *))lcg64_default;
	random->write  = (char* (*)(void *))lcg64_write;
	random->bread  = (unsigned (*)(unsigned char *,void *))lcg64_bread;
	random->bwrite  =(unsigned (*)(unsigned char *,void *))lcg64_bwrite;
	random->uniformRandom  =(double (*)(void *))lcg64;
	random->uniformRandom2  =(TWO_VECTOR (*)(void *))lcg64_2;
}
MPI_Datatype lcg64_MPIType(void)
{
	static MPI_Datatype PARM_TYPE;
	static int initialized = 0;
	if (initialized == 0)
	{
		MPI_Aint disp[3]; 
		int blkcnt[] = { 1, 1, 1, };
		MPI_Datatype types[] = {MPI_LONG_LONG, MPI_INT, MPI_INT,  };
		int n = 3;
		LCG64_PARM parm;
        MPI_Get_address(&(parm.state), disp + 0);
        MPI_Get_address(&(parm.multID), disp + 1);
        MPI_Get_address(&(parm.prime), disp + 2);
		for (int i = 1; i < n; i++) disp[i] -= disp[0];
		disp[0] = 0;
        MPI_Type_create_struct(n, blkcnt, disp, types, &PARM_TYPE);
		MPI_Type_commit(&PARM_TYPE);
		initialized = 1;
	}
	return PARM_TYPE; 
}

char *lcg64_write(LCG64_PARM *parm)
{
	static char line[64];
   if (parm != NULL) assert(parm->multID < 3); 
	if (parm != NULL)  sprintf(line, "%16.16"PRIx64" %1u %8.8x", parm->state, parm->multID, parm->prime);
	else               sprintf(line, "%16.16llx %1u %8.8x", 0LL,0,0);
	return line;
}

unsigned lcg64_bread(unsigned char* in, LCG64_PARM *parm)
{
   if (in != NULL)
   {
      parm->state = mkInt(in, "u8");
      parm->multID = mkInt(in+8, "u4"); 
      parm->prime = mkInt(in+12, "u4");
   }
   return 16;
}
unsigned lcg64_bwrite(unsigned char* out, LCG64_PARM *parm)
{
   if (out != NULL)
   {
      copyBytes(out,    &parm->state,  8);
      copyBytes(out+8,  &parm->multID, 4);
      copyBytes(out+12, &parm->prime,  4);
   }
   return 16;
}

void lcg64_write_stdout(LCG64_PARM *parm)
{
	printf("write_stdout %d: %"PRIx64" %x %x\n", getRank(0), parm->state, parm->multID, parm->prime);
}

char *lcg64_parse(char *string, LCG64_PARM *parm)
{
	int cnt,charCnt; 
	cnt=sscanf(string, "%"PRIx64" %u %x%n", &parm->state, &parm->multID, &parm->prime,&charCnt);
	if (cnt != 3 || parm->prime==0 || parm->multID > 2) charCnt=0; 
	return string+charCnt;
}

void lcg64_default(unsigned n, LONG64 seed, LCG64_PARM *parm)
{
	static int multID = 0;
	static LONG64 prime = 0ll;
	if (parm == NULL) return; 

	if (multID == 0) prime = nextPrime(); 

	parm->multID = multID;
	parm->state = INIT_SEED ^ seed ;
	parm->prime = prime;
	multID = (multID+1 )%3; 
}
int lcg64_checkValue(LCG64_PARM *parm)
{
	unsigned int multID = parm->multID ;
	LONG64 state = parm->state ;
	LONG64 prime = parm->prime ;
   if (multID > 2) return 1; 
   if (state == 0) return 1; 
   if (prime%2 == 0 ) return 1; 
   return 0; 
}

double lcg64(LCG64_PARM *parm)
{
	double r;
	LONG64 MULT[] = { 0x27bb2ee687b0b0fdllu, 0x2c6fe96ee78b6955llu, 0x369dea0f31a53f85llu };
	parm->state = MULT[parm->multID]*parm->state + parm->prime;
	r = parm->state*TWO_M64;
	return r;
}
TWO_VECTOR lcg64_2(LCG64_PARM *parm)
{
	TWO_VECTOR r;
	LONG64 MULT[] = { 0x27bb2ee687b0b0fdllu, 0x2c6fe96ee78b6955llu, 0x369dea0f31a53f85llu };
	parm->state = MULT[parm->multID]*parm->state + parm->prime;
	r.x = parm->state*TWO_M64;
	parm->state = MULT[parm->multID]*parm->state + parm->prime;
	r.y = parm->state*TWO_M64;
	return r;
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
