/* Version 1 of Reflect                                    */
/*                                                         */
/*  hard reflection from top and bottom of box   `          */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "system.h"
#include "ddcMalloc.h"
typedef struct reflect_parms_str
{
	int dummy;
} REFLECT_PARMS;

void fsum(double *, double *, double *, THREE_MATRIX *);
int getRank(int);
int nelements;

/* REFLECT potential routines */
REFLECT_PARMS *reflect_parms(OBJECT*object)
{
	REFLECT_PARMS *parms;
	parms = ddcMalloc(sizeof(REFLECT_PARMS));
	return parms;
}

RCUT_TYPE *reflectCutoff(SYSTEM* sys, REFLECT_PARMS*parms, int *n)
{
	static RCUT_TYPE rcut[2];
	*n = 1;
	rcut[0].value = 0.0;
	rcut[0].mode = RCUT_ALL;
	rcut[1].value = 0.0;
	return rcut;
}


void reflect(SYSTEM*sys, REFLECT_PARMS*parms, ETYPE*e)
{
	int nlocal;
	double *vz, *rz;
	double topedge,botedge;
	void box_get(BOX_STRUCT *, int, void *);
	THREE_MATRIX *h0;
	THREE_VECTOR *c0;
	/* REFLECT vars */
	/* Initialize Variables  */

	nlocal = sys->nlocal;
	box_get(NULL, HO_PTR, &h0);
	box_get(NULL, CORNER_PTR, &c0);

        rz = sys->collection->state->rz; 
        vz = sys->collection->state->vz; 
	botedge=c0->z;
	topedge=botedge+h0->zz;
	for (int i = 0; i < nlocal; i++)
	{
		if(rz[i]>topedge)
		{
			vz[i]=-vz[i];
			rz[i]= 2.0*topedge-rz[i];
		}
		else if(rz[i]<botedge)
		{
			vz[i]=-vz[i];
			rz[i]= 2.0*botedge-rz[i];
		}
	}

	
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
