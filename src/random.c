#include "random.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "three_algebra.h"
#include "object.h"
#include "ddcMalloc.h"
#include "lcg64.h"
#include "pio.h"
#include "particle.h"
#include "utilities.h"
#include "mpiUtils.h"
#include "object.h"
#include "kn.h"

/**  Generates a 64-bit random seed based on the time of day and rank in
 *   MPI_COMM_WORLD.
 *
 *   Notes: rand_r actually returns only 31 random bits so bits 32 and
 *   64 are not actually random in this implementation.  Oh well.
 */
LONG64 generateRandomSeed(void)
{
   LONG64 epoch = getTick();
   unsigned seed1 = epoch & 0x55555555;
   unsigned seed2 = epoch & 0xaaaaaaaa;
   epoch = epoch >> 32;
   seed1 |= epoch & 0xaaaaaaaa;
   seed2 |= epoch & 0x55555555;

   seed1 ^= getRank(0);
   seed2 ^= getRank(0);
   for (unsigned ii=0; ii<1000; ++ii)
   {
      rand_r(&seed1);
      rand_r(&seed2);
   }
   LONG64 seed = rand_r(&seed2);
   seed = (seed << 32) + rand_r(&seed1);
   return seed; 
}

RANDOM *random_init(void *parent, char *name)
{
	static RANDOM *random;
	int randomizeSeed; 
	char *type;
	random = (RANDOM *) object_initialize(name, "RANDOM", sizeof(RANDOM));
	object_get((OBJECT *) random, "type", &(type), STRING, 1, NULL);
	object_get((OBJECT *) random, "seed", &random->seed, U64, 1, "0");
	object_get((OBJECT *) random, "randomizeSeed", &randomizeSeed, INT, 1, "0");
	if (randomizeSeed == 1) random->seed = generateRandomSeed(); 
	random->parent = parent; 
	random->name = strdup(name);
	random->type = strdup(type);
	random->parmsArray = NULL;
   random->useDefault=0; 
	if (strcmp(type,"LCG64")==0) 
   {
      lcg64_init(random);
   }
   random_registerParms(random);

	return random;
}
void random_registerParms(RANDOM* random)
{
	static int registered = 0;
	if (registered != 0) return;

	MPI_Datatype PARM_TYPE=(MPI_Datatype)0LL; 
	if (random->itype==LCG64)
		PARM_TYPE = lcg64_MPIType(); 
	particleRegisterinfo((void *)&(random->parmsArray), random->sizeOfParms, random->sizeOfParms, PARM_TYPE, NULL, NULL);
	registered = 1;
}
void *random_getParms(RANDOM *random, unsigned k)
{
   if (random == NULL) return NULL; 
	return (char *)(random->parmsArray)+(random->sizeOfParms*k); 
}
int random_checkValue(RANDOM *random, unsigned k)
{
   int flag=0; 
   void *value = (random->parmsArray)+(random->sizeOfParms*k); 
	if (random->itype==LCG64) flag = lcg64_checkValue(value);
	return flag; 
}

/** Creates a new default randomParms including allocating memory */
void* random_newParms(RANDOM* this)
{
	void* parms = ddcMalloc(this->sizeOfParms);
	this->defaultValue(0, this->seed, parms);
	return parms;
}

double gasdev0(PRAND48_STATE* handle)
{
	double fac, rsq, v1, v2;
	do
	{
		v1 = (2.0*prand48(handle) - 1.0);
		v2 = (2.0*prand48(handle) - 1.0);
		rsq = v1*v1 + v2*v2;
	}
	while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	return v2*fac;
}
THREE_VECTOR gasdev3d0(double sigma, PRAND48_STATE* handle)
{
	THREE_VECTOR r;
	r.x = sigma*gasdev0(handle); 
	r.y = sigma*gasdev0(handle); 
	r.z = sigma*gasdev0(handle); 
	return r; 
}

double gasdev(RANDOM *random,void *parm)
{
	double fac, rsq, v1, v2;
	do
	{
		v1 = (2.0*random->uniformRandom(parm) - 1.0);
		v2 = (2.0*random->uniformRandom(parm) - 1.0);
		rsq = v1*v1 + v2*v2;
	}
	while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	return v2*fac;
}
THREE_VECTOR gasdev3d(RANDOM *random,void *parm)
{
	THREE_VECTOR r; 
   TWO_VECTOR v; 
	double fac, rsq ;
	do
	{
		v = random->uniformRandom2(parm);
		v.x = (2.0*v.x - 1.0);
		v.y = (2.0*v.y - 1.0);
		rsq = v.x*v.x + v.y*v.y;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	r.x = v.x*fac; 
	r.y = v.y*fac; 
	do
	{
		v = random->uniformRandom2(parm);
		v.x = (2.0*v.x - 1.0);
		v.y = (2.0*v.y - 1.0);
		rsq = v.x*v.x + v.y*v.y;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	r.z = v.x*fac;
	return r; 
}


/** Relativistic Maxwell-Juttner random number generator
 *  by Ian Ellis.
 *  NOTE: T is in terms of the rest energy of the species
 *  (i.e. send kB*T/(m*c^2)) and this returns p = gamma*beta
 */
THREE_VECTOR jutdev0(PRAND48_STATE* handle, double T)
{
    static double TBessK=0.0, jutPeak=0.0, jutMax=0.0, a0=0.0, a0sq=0.0;
    double rand1, rand2, rand3, rsq, rr;
    double gamma, magP, rLorentzVal, jutVal;
    THREE_VECTOR p;

    //misc. stuff we only need to do once during the simulation
    if (TBessK==0.0)
    {
        TBessK = T*kn(2,1.0/T);
        if (TBessK==0.0)
        {
            printf("FATAL: The Juttner distribution generator\n"
            "       cannot support such low temperatures.  Please\n"
            "       use a higher temperature or a Maxwellian distribution.\n\n");
            abortAll(-1);
        }

        //the location of the peak of the PDF
        jutPeak = creal(cpow((16.0*T*T*T+3.0*sqrt(3.0)*csqrt(-32.0*pow(T,4)-13.0*T*T-4.0)-9.0*T),1.0/3.0)/(3.0*pow(2.0,1.0/3.0))
            +(pow(2.0,1.0/3.0)*(4.0*T*T+3.0))/(3.0*cpow(16.0*T*T*T+3.0*sqrt(3.0)*csqrt(-32.0*pow(T,4)-13.0*T*T-4.0)-9.0*T,1.0/3.0))
            +(2.0*T)/3.0);

        //the maximum of the PDF
        jutMax = jutPeak*sqrt(jutPeak*jutPeak-1)*exp(-jutPeak/T)/TBessK;

        a0 = 2.0*T;
        a0sq = a0*a0;
    }

    //find gamma using the rejection method, testing against a Lorentzian
    do
    {
        do
        {
            do
            {
                rand1 = prand48(handle);
                rand2 = 2.0*prand48(handle)-1.0;
            }   // for some reason catching the 0.0 isn't SOP...
            while (rand1*rand1 + rand2*rand2>1.0 || rand1==0.0);

            gamma = a0*rand2/rand1 + jutPeak;
        }
        while (gamma<1.0);

        rLorentzVal = (1+(gamma-jutPeak)*(gamma-jutPeak)/a0sq)/(1.1*jutMax);
        jutVal = gamma*sqrt(gamma*gamma-1.0)*exp(-gamma/T)/TBessK;
    }
    while (prand48(handle) > jutVal*rLorentzVal);

    //determine the directon this thing points
    do
    {
        rand1 = 2.0*prand48(handle)-1.0;
        rand2 = 2.0*prand48(handle)-1.0;
        rand3 = 2.0*prand48(handle)-1.0;
        rsq = rand1*rand1+rand2*rand2+rand3*rand3;
    }
    while (rsq>1.0 || rsq==0.0);

    //here's the momentum, finally!
    magP = sqrt((gamma*gamma)-1);

    rr = 1/sqrt(rsq);
    p.x = rand1*magP*rr;
    p.y = rand2*magP*rr;
    p.z = rand3*magP*rr;

    return p;
}


THREE_VECTOR jutdev(RANDOM *random, void *parm, double T)
{
    static double TBessK=0.0, jutPeak=0.0, jutMax=0.0, a0=0.0, a0sq=0.0;
    double rand1, rand2, rand3, rsq, rr;
    double gamma, magP, rLorentzVal, jutVal;
    THREE_VECTOR p;

    //misc. stuff we only need to do once during the simulation
    if (TBessK==0.0)
    {
        TBessK = T*kn(2,1.0/T);
        if (TBessK==0.0)
        {
            printf("FATAL: The Juttner distribution generator\n"
            "       cannot support such low temperatures.  Please\n"
            "       use a higher temperature or a Maxwellian distribution.\n\n");
            abortAll(-1);
        }

        //the location of the peak of the PDF
        jutPeak = creal(cpow((16.0*T*T*T+3.0*sqrt(3.0)*csqrt(-32.0*pow(T,4)-13.0*T*T-4.0)-9.0*T),1.0/3.0)/(3.0*pow(2.0,1.0/3.0))
            +(pow(2.0,1.0/3.0)*(4.0*T*T+3.0))/(3.0*cpow(16.0*T*T*T+3.0*sqrt(3.0)*csqrt(-32.0*pow(T,4)-13.0*T*T-4.0)-9.0*T,1.0/3.0))
            +(2.0*T)/3.0);

        //the maximum of the PDF
        jutMax = jutPeak*sqrt(jutPeak*jutPeak-1)*exp(-jutPeak/T)/TBessK;

        a0 = 2.0*T;
        a0sq = a0*a0;
    }

    //find gamma using the rejection method, testing against a Lorentzian
    do
    {
        do
        {
            do
            {
                rand1 = random->uniformRandom(parm);
                rand2 = 2.0*random->uniformRandom(parm)-1.0;
            } // for some reason catching the 0.0 isn't SOP...
            while (rand1*rand1 + rand2*rand2>1.0 || rand1==0.0);

            gamma = a0*rand2/rand1 + jutPeak;
        }
        while (gamma<1.0);

        rLorentzVal = (1+(gamma-jutPeak)*(gamma-jutPeak)/a0sq)/(1.1*jutMax);
        jutVal = gamma*sqrt(gamma*gamma-1.0)*exp(-gamma/T)/TBessK;
    }
    while (random->uniformRandom(parm) > jutVal*rLorentzVal);

    //determine the directon this thing points
    do
    {
        rand1 = 2.0*random->uniformRandom(parm)-1.0;
        rand2 = 2.0*random->uniformRandom(parm)-1.0;
        rand3 = 2.0*random->uniformRandom(parm)-1.0;
        rsq = rand1*rand1+rand2*rand2+rand3*rand3;
    }
    while (rsq>1.0 || rsq==0.0);

    //here's the momentum, finally!
    magP = sqrt((gamma*gamma)-1);

    rr = 1/sqrt(rsq);
    p.x = rand1*magP*rr;
    p.y = rand2*magP*rr;
    p.z = rand3*magP*rr;

    return p;
}


/** This is hardwired to create a LCG64 type generator.  If we ever
 *  implement another type, and if the clients of disposable generators
 *  care, we can add a type parameter.
 */
void random_getDisposableRandomGenerator(
	LONG64 seed, RANDOM** random, void** parms)
{
	*random = ddcMalloc(sizeof(RANDOM));
	(*random)->name = NULL;
	(*random)->objclass = NULL;
	(*random)->value = NULL;
	(*random)->type = NULL;
	(*random)->parent = NULL;
	(*random)->parmsArray = NULL;
	(*random)->seed = seed;
	lcg64_init(*random);
	*parms = random_newParms(*random);
}

/** This function is intended to be called only for disposable
 *  generators.  For disposables none of the pointers in the RANDOM
 *  structure were allocated so there is no need to free them.
 */
void random_free(RANDOM* this, void* parms)
{
	ddcFree(parms);
	ddcFree(this);
}

/** Random number generator that combines a particle gid, a seed, and a
 *  call site id to generate a random number.  Since the underlying
 *  random engine is erand48, only the 48 low order bits in the three
 *  parameters are used.
 *
 *  The rationale for generating random numbers this way is that for a
 *  fixed user specified seed, we want reproducible results, regardless
 *  of the number of tasks the code is running on.  Since we generate a
 *  new seed for erand based on the user seed and the particle gid, we
 *  get the same number (for a given seed) regardless of the calling
 *  history of erand.  The callSite paremeter allows for different call
 *  sites in the code to get different pseudo-random streams even for
 *  the same particle and user seed.
 *
 *  Issue: To my surprise, the random streams for nearly equal gids tend
 *  to be correlated, at least for a while.  This doesn't hurt anything
 *  when the gids are randomly distributed to the particles, but could
 *  have a nasty effect if there are spatial correlations between the
 *  gids.
*/
double prand48(PRAND48_STATE* this)
{
   return erand48(this->seedVec);
}

PRAND48_STATE prand48_init(LONG64 gid, LONG64 seed, LONG64 callSite)
{
	PRAND48_STATE handle;
	
   seed ^= callSite;
   gid ^= seed;
	
   handle.seedVec[0] = handle.seedVec[1] = handle.seedVec[2] = 0;
   unsigned short mask[4] = {0x000F, 0x00F0, 0x0F00, 0xF000};
   for (unsigned ii=0; ii<4; ++ii)
   {
      handle.seedVec[0] += (gid & mask[ii]); gid = gid >> 4;
      handle.seedVec[1] += (gid & mask[ii]); gid = gid >> 4;
      handle.seedVec[2] += (gid & mask[ii]); 
   }
	return handle;
}
void testForRandomKeyword(OBJECT *obj)
{
   if (object_testforkeyword(obj, "random"))
   {
      if (getRank(0) == 0) 
      {
         printf("Keyword random now needs to be define in the system object.  Remove keyword %s from %s object and added to SYSTEM object\n","random", obj->name); 
         exit(1); 
      }
   }
}

void missingRandomError(char *string)
{
   if (getRank(0) != 0) return;
   printf("ERROR: A  %s has failed to find a valid random OBJECT\n"
         "       You probably forgot to specify the random keyword in SYSTEM object\n",string);
   abortAll(-1);
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
