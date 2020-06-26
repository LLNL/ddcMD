#ifndef PARTICLE_H
#define PARTICLE_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef  void (*initType)(void *, void *start, int size, int stride, int begin, int n) ;
typedef struct particleinfo_st
{
	void **start;
	int size, stride,nParticle;
	MPI_Datatype type;
   initType init; 
   void *initParms; 
} PARTICLEINFO;

void particleRegisterinfo(void **start, int size, int stride, MPI_Datatype type, initType init, void *initParms);
void particleGetinfo(void *p, int element);
void particleSort(int *map, int n);
void particleRegisterAccum(void **start, int size, int stride, MPI_Datatype type, initType, void *initParms);
void particlePrintinfo(char *label);
void particleAccum(double *recv, int nlist, int *list);
void particleFillBuffer(void *buffer, int first, int n);
void particlePutinfo(void *ptr, int element);
int particleSizeAccum(void);
int particleSizeinfo(void);
int particleAllocateinfo(int nparticles);
void particleSortinfo(char *map, int mapstride, int n);/*Sorts particles by map */ 
#ifdef __cplusplus
}
#endif

#endif // #ifndef PARTICLE_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
