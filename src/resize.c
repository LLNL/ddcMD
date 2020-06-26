#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "three_algebra.h"
#include "collection.h"
#include "ddcMalloc.h"
#include "pio.h"

static int maxsize1 = 0;
static int maxsize2 = 0;
static int infosize;
static int size = 0;
void mgpt2ddcmemset(STATE*state);
int particleAllocateinfo(int n);
void resizePrint(PFILE*file)
{
	Pprintf(file, "resize controls %d kBytes of memory\n", size/1024);
	//fflush(file);
}

void resize(int n, int mode, STATE*state)
{
  const int nmin = 5000/*50000*/;
	if (n < nmin) n = nmin;
	switch (mode)
	{
	case 1:
		if (n < maxsize1) return;
		maxsize1 = n;
		state->rx = ddcRealloc(state->rx, n*sizeof(double));
		state->ry = ddcRealloc(state->ry, n*sizeof(double));
		state->rz = ddcRealloc(state->rz, n*sizeof(double));
		state->label = ddcRealloc(state->label, n*sizeof(gid_type));
		state->atomtype = ddcRealloc(state->atomtype, n*sizeof(int));
		state->species = ddcRealloc(state->species, n*sizeof(SPECIES *));
		state->group = ddcRealloc(state->group, n*sizeof(GROUP *));
		particleAllocateinfo(n);
		break;
	case 2:
		if (n < maxsize2) return;
		maxsize2 = n;
		state->rx = ddcRealloc(state->rx, n*sizeof(double));
		state->ry = ddcRealloc(state->ry, n*sizeof(double));
		state->rz = ddcRealloc(state->rz, n*sizeof(double));
		state->label = ddcRealloc(state->label, n*sizeof(gid_type));
		state->atomtype = ddcRealloc(state->atomtype, n*sizeof(int));
		state->species = ddcRealloc(state->species, n*sizeof(SPECIES *));
		state->group = ddcRealloc(state->group, n*sizeof(GROUP *));
		state->fx = ddcRealloc(state->fx, n*sizeof(double));
		state->fy = ddcRealloc(state->fy, n*sizeof(double));
		state->fz = ddcRealloc(state->fz, n*sizeof(double));
		state->q = ddcRealloc(state->q, n*sizeof(double));
		state->potentialEnergy = ddcRealloc(state->potentialEnergy, n*sizeof(double));
		state->virial = ddcRealloc(state->virial, n*sizeof(double));
		state->sion = ddcRealloc(state->sion, n*sizeof(THREE_SMATRIX));
		infosize = particleAllocateinfo(n);
		size = infosize + n*(sizeof(double)*8 + sizeof(int)*1 + sizeof(gid_type) + sizeof(SPECIES *)*2);
		break;
      
	}
	mgpt2ddcmemset(state);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
