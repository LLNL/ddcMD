#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include "three_algebra.h"
//#include "ddc.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#define MAXSIZE  512
#define MAXELEMENTS 128

/*static PARTICLEINFO *info = NULL;*/
static PARTICLEINFO info[MAXELEMENTS];
static int nelements = 0;
static int particlesize = 0;

//static PARTICLEINFO *accum = NULL;
static PARTICLEINFO accum[MAXELEMENTS];
static int naccum = 0;
static int accumsize = 0;

void *ExpandBuffers(void *, int, int, int, char *, char *);

void particleRegisterinfo(void **start, int size, int stride, MPI_Datatype type, initType init, void *initParms)
{
	int mpisize;
/*	info = ExpandBuffers(info, sizeof(PARTICLEINFO), nelements + 1, 8, LOCATION("particleReisterinfo"),"info");*/
	MPI_Type_size(type, &mpisize);
	info[nelements].start = start;
	info[nelements].size = size;
	info[nelements].stride = stride;
	info[nelements].type = type;
	info[nelements].init = init;
	info[nelements].initParms = initParms;
	info[nelements].nParticle = 0;
	particlesize += size;
	nelements++;
	if (nelements > MAXELEMENTS)
	{
		printf("Exceed MAX number of elements in particleRegisterinfo\n");
		MPI_Finalize();
	}
	if (particlesize > MAXSIZE)
	{
		printf("Exceed MAX size in particleRegisterinfo\n");
		MPI_Finalize();
	}
}

void particleRegisterAccum(void **start, int size, int stride, MPI_Datatype type, initType init, void *initParms)
{
	int mpisize;
	MPI_Type_size(type, &mpisize);
	info[nelements].start = start;
	if (type != MPI_DOUBLE)
	{
		printf("Error in particleRegisterAccum\n");
		abortAll(1);
	}
	//accum = ExpandBuffers(accum, sizeof(PARTICLEINFO), nelements + 1, 8,LOCATION("particleRegisterAccum"), "accum");
	accum[naccum].start = start;
	accum[naccum].size = size;
	accum[naccum].stride = stride;
	accum[naccum].type = type;
	accumsize += size;
	naccum++;
	if (naccum > MAXELEMENTS)
	{
		printf("Exceed MAX number of elements in particleRegisteraccum\n");
		abortAll(1);
	}
	if (accumsize > MAXSIZE)
	{
		printf("Exceed MAX size in particleRegisteraccum\n");
		abortAll(1);
	}
}

void particlePrintinfo(char *label)
{
	void **start;
	int i, size;
	for (i = 0; i < nelements; i++)
	{
		start = info[i].start;
		size = info[i].size;
		printf("%s i=%d start=%p size=%d\n", label, i, start, size);
		fflush(stdout);
	}
}

int particleAllocateinfo(int nParticle)
{
	void *start;
	int i, size, stride, cnt;
	cnt = 0;
   
   int nAllocate = nParticle;
   if (nAllocate < 5000) nAllocate =5000;
	for (i = 0; i < nelements; i++)
	{
		start = *info[i].start;
      int begin = info[i].nParticle;
		size = info[i].size;
		stride = info[i].stride;
		start = ddcRealloc(start, size*nAllocate);
		*info[i].start = start;
		cnt += size*nAllocate;
      if (nParticle > info[i].nParticle)
      {
            if (info[i].init != NULL) info[i].init(info[i].initParms,start,size,stride,begin,nParticle); 
      }
      info[i].nParticle = nParticle; 
	}
	return cnt;
}

void *GetLocationinfo(int n, int index, MPI_Datatype*TYPE)
{
	MPI_Aint displs[MAXELEMENTS];
	MPI_Datatype htype[MAXELEMENTS];
	int i, blockcount[MAXELEMENTS];
	char *hstart;
	assert(n>0);
	for (i = 0; i < nelements; i++)
	{
		MPI_Datatype type = info[i].type;
		MPI_Aint stride = info[i].stride;
		MPI_Type_create_hvector(n, 1, stride, type, htype + i);
		MPI_Type_commit(htype+i); 
		hstart = *info[i].start;
		hstart += stride*index;
		MPI_Get_address(hstart, &displs[i]);
		blockcount[i] = 1;
	}

	for (i = nelements; i > 0; i--) 
	  displs[i - 1] -= displs[0];

	if (*(TYPE) != 0) MPI_Type_free(TYPE);
	MPI_Type_create_struct(nelements, blockcount, displs, htype, TYPE);
	MPI_Type_commit(TYPE);
	hstart = *info[0].start;
	hstart += info[0].stride*index;
	for (i = 0; i < nelements; i++) MPI_Type_free(htype+i); 
	return hstart;
}

void *GetLocationAccum(int n, int index, MPI_Datatype*TYPE)
{
	MPI_Aint displs[MAXELEMENTS];
	MPI_Datatype type, htype[MAXELEMENTS];
	int i, blockcount[MAXELEMENTS];
	MPI_Aint stride;
	char *hstart;
	for (i = 0; i < naccum; i++)
	{
		type = accum[i].type;
		stride = accum[i].stride;
		assert(n>0);
		MPI_Type_create_hvector(n, 1, stride, type, htype + i);
		hstart = *accum[i].start;
		hstart += stride*index;
		MPI_Get_address(hstart, &displs[i]);
		blockcount[i] = 1;
	}

	for (i = naccum; i > 0; i--)
		displs[i - 1] -= displs[0];
	if (*(TYPE) != 0) MPI_Type_free(TYPE);
	MPI_Type_create_struct(naccum, blockcount, displs, htype, TYPE);
	MPI_Type_commit(TYPE);
	hstart = *accum[0].start;
	hstart += accum[0].stride*index;
	return hstart;
}

int particleSizeinfo(void)
{
	return particlesize;
}

int particleSizeAccum(void)
{
	return accumsize;
}

void particleAccum(double *recv, int nlist, int *list)
{
	char *ptr;
	for (int i = 0; i < naccum; i++)
	{
		int stride = accum[i].stride;
		//int size = accum[i].size;
		ptr = *accum[i].start;
		for (int j = 0; j < nlist; j++) *(double *)(ptr + stride*list[j]) += recv[j];
		recv += nlist;
	}
}

void particleFillBuffer(void *buffer, int first, int n)
{
	int i, j, stride, size;
	char *ptrA, *ptrB;
	ptrA = buffer;
	for (i = 0; i < nelements; i++)
	{
		stride = info[i].stride;
		size = info[i].size;
		ptrB = *info[i].start;
		ptrB += stride*first;
		for (j = 0; j < n; j++)
		{
			memcpy(ptrA, ptrB, size);
			ptrA += size;
			ptrB += stride;
		}
	}
}

void particleGetinfo(void *p, int element) /*Gets particle information */
{
	int i, k;
	char *start;
	k = 0;
	for (i = 0; i < nelements; i++)
	{
		start = *info[i].start;
		memcpy(((char *)p) + k, start + info[i].stride*element, info[i].size);
		k += info[i].size;
	}
}

void particlePutinfo(void *ptr, int element)
{
	int i;
	char *start;
	char *p;
	p = ptr;
	for (i = 0; i < nelements; i++)
	{
		start = *info[i].start;
		memcpy(start + info[i].stride*element, (void *)p, info[i].size);
		p += info[i].size;
	}
}

void particleSortinfo(char *map, int mapstride, int n) /*Sorts particles by map */ 
{
	int i, j, js, k;
	char q[MAXSIZE], s[MAXSIZE];
	js = 0;
	for (k = 0; k < n;)
	{
		for (j = js; j < n; j++) if (*(int *)(map + mapstride*j) != j) break;
		if (j == n) break;
		js = j;
		i = *(int *)(map + mapstride*j);
		particleGetinfo(s, js);
		while (i != js)
		{
			particleGetinfo(q, i);
			particlePutinfo(q, j);
			*(int *)(map + mapstride*j) = j;
			j = i;
			i = *(int *)(map + mapstride*j);
			k++;
		}
		k++;
		*(int *)(map + mapstride*j) = j;
		particlePutinfo(s, j);
		js++;
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
