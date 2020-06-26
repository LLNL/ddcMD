#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(BGL) || defined(BGP) 
#include <mass.h> 
#include <massv.h> 
#endif 
#include <string.h>
#include <mpi.h>
#include "three_algebra.h"
#include "species.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "ptiming.h"
#include "opt.h"
#include "species.h"
#include "system.h"
#include "box.h"
#include "particle.h"
#include "opt.h"
#include "ddcMalloc.h"
#include "heap.h"
#include "units.h"
#include "hpmWrapper.h"

typedef struct PAIRENERGY_parms_str
{
	double rmax,r_expansion;
	int nspecies; 
	void *sparms; 
	int cmax; 
} PAIRENERGY_PARMS;

typedef  struct pair_pass_parms_st { double v[72],d[72];} PAIR_PASS_PARMS ;
int np,nspecies; 
double rcut,r2cut; 
double	sion_xx,sion_xy,sion_xz,sion_yy,sion_yz,sion_zz ;
int offset[4]; 
extern FILE *ddcfile;
double  pairKernel		(void *parms, OPT_BOX *box, OPT_PARTICLE *rv, THREE_VECTOR *s, unsigned *type, double *dF, THREE_VECTOR *f);

static double  (*kernel)    		(void *parms, OPT_BOX *box, OPT_PARTICLE *rv, THREE_VECTOR *s, unsigned *type, double  *dF, THREE_VECTOR *f);




static  OPT_PARTICLE *rv=NULL; 
static  double *dF=NULL; 
static OPT_BOX *boxlist; 
static int nbox; 
void pairEnergy_start(PAIRENERGY_PARMS *parms) 
{
	np=0; 
	nspecies=parms->nspecies; 
	rcut = parms->rmax;  
	r2cut = rcut*rcut; 
}
void pairfs_parms(POTENTIAL *object, PAIRENERGY_PARMS *parms)
{
	int j,nspecies;
	PAIR_PASS_PARMS **pp; 
	char **namelist,*mode; 

	nspecies = parms->nspecies; 
	species_get(NULL,NAMELIST,(void *)&namelist);
//	pp = (PAIR_PASS_PARMS **)parms->pass_parms; 
	pp = NULL; 
	for (int i=0;i<nspecies;i++) 
	{
		if (pp[i+i*nspecies]==NULL)
		   ddcMallocAligned((void*)&pp[i+i*nspecies],16,sizeof(PAIR_PASS_PARMS));
		for (j=i+1;j<nspecies;j++) 
		{
			if (pp[i+j*nspecies]==NULL)
			   ddcMallocAligned((void*)&pp[i+j*nspecies],16,sizeof(PAIR_PASS_PARMS));
			pp[j+i*nspecies ] = pp[i+j*nspecies];
		}
	}
	object_get((OBJECT *) object, "mode", &mode, STRING, 1, "ANALYTIC"); 
}

PAIRENERGY_PARMS *pairEnergy_parms(POTENTIAL *object)
{
	PAIRENERGY_PARMS *parms;
	int nspecies;
	char **namelist; 
	species_get(NULL,NAMELIST,(void *)&namelist);
	nspecies=system_getNspecies((SYSTEM *)object->parent);

	parms = malloc(sizeof(PAIRENERGY_PARMS));
	parms->nspecies=nspecies; 
//	parms->pass_parms=calloc(nspecies*nspecies,sizeof(void*)); 

//	object_get((OBJECT *) object, "rmax", &parms->rmax, WITH_UNITS, 1, "0.0","Angstrom",NULL);
// object_get((OBJECT *) object, "r_expansion", &parms->r_expansion, WITH_UNITS, 1, "3.0","Angstrom",NULL); 
	pairfs_parms(object,parms);

	kernel = pairKernel ;   
	return parms;
}
RCUT_TYPE *pairEnergyCutoff(SYSTEM*sys, PAIRENERGY_PARMS*parms, int *n)
{
	static RCUT_TYPE rcut[2];
	static int ncut = 1;
	rcut[0].value = parms->rmax;  
	rcut[0].mode = RCUT_ALL;
	rcut[ncut].value = parms->rmax;
	rcut[ncut].mode = RCUT_ALL;
	*n = ncut;
	return rcut;
}
void pairEnergy(SYSTEM *sys,PAIRENERGY_PARMS *parms,ETYPE *energyInfo )
{
	static unsigned *sort_index=NULL, *type=NULL; 
	static THREE_VECTOR *f=NULL; 
	OPT_BOX_SIZE box_size; 
	THREE_VECTOR dia,shift_zero = {0.0,0.0,0.0},shift_minus,shift_plus; 
	int i,j,nion,nlocal;
	int ix,iy,iz,nx,ny,nz; 
	double *rx,*ry,*rz,*fx,*fy,*fz; 
	double rmax,pair_energy; 
	int nion_16;
	unsigned int block_number; 
	
	HPM_Start("pairEnergy");
	char* scratch = (char *)heapGet(&block_number); 
	rmax = parms->rmax;
 	dia=box_get_diagonal(NULL);
	nion = sys->nion; 
	nlocal = sys->nlocal; 
	rx = sys->collection->state->rx; 
	ry = sys->collection->state->ry; 
	rz = sys->collection->state->rz; 
	fx = sys->collection->state->fx; 
	fy = sys->collection->state->fy; 
	fz = sys->collection->state->fz; 

	nion_16= (nion+128)&0xfffffff0;
	size_t size = 0;
	rv         = (OPT_PARTICLE *)(scratch+size);  size += sizeof(OPT_PARTICLE) * 2*nion_16;
	sort_index = (unsigned *)(scratch+size);      size += sizeof(int)*nion_16; 
	type 	     = (unsigned *)(scratch+size);      size += sizeof(unsigned)*nion_16;
	f          = (THREE_VECTOR *)(scratch+size);  size += sizeof(THREE_VECTOR)*nion_16; 
	dF         = (double *)(scratch+size);        size += sizeof(double )*nion_16; 

	heapEndBlock(block_number, size);
	

	HPM_Start("opt_trans");
	nion=opt_transform(rcut,nion,nlocal,rx, ry, rz, rv) ;
	HPM_Stop("opt_trans");

	HPM_Start("opt_box");
	boxlist = opt_box(rmax, nion, nlocal, rv, sort_index, &box_size);
	nbox = box_size.nbox; 
	nx = box_size.nx; 
	ny = box_size.ny; 
	nz = box_size.nz; 
	offset[0] = nx -1; 
	offset[1] = (ny-1)*nx -1; 
	offset[2] = (ny)*nx -1; 
	offset[3] = (ny+1)*nx -1; 
	HPM_Stop("opt_box");

	shift_zero.x = shift_zero.y = shift_zero.z = 0.0; 
	shift_plus.x = dia.x ;    shift_plus.y =0.0; shift_plus.z  =0.0; 
	shift_minus.x = -dia.x ; shift_minus.y =0.0; shift_minus.z =0.0; 
	pairEnergy_start(parms); 
	pair_energy=sion_xx=sion_xy=sion_xz=sion_yy=sion_yz=sion_zz=0.0 ;

	for (i=0;i<nion;i++) 
	{
		j = sort_index[i]; 
		if (j < nlocal )  dF[i] =0.5 ;  else  dF[i] =0.0; 
		f[i].x=f[i].y=f[i].z=0.0; 
		type[i] = sys->collection->state->atomtype[j] & 0xffff; 
	}
	HPM_Start("block0");
	for (iz=0;iz<nz;iz++) 
	{
		i =  1 + nx*(1 + ny *iz); 
		if (boxlist[i].edge == 0 ) 
		{
			for (iy=1;iy<ny-1;iy++) 
			{
				i = 1 + nx*(iy + ny *iz); 
				pair_energy += kernel(parms, boxlist+i, rv, &shift_plus, type,  dF,f);i++;
			for (ix=2;ix<nx-2;ix++) {pair_energy += kernel(parms, boxlist+i, rv, &shift_zero, type,  dF,f);i++;}
				pair_energy += kernel(parms, boxlist+i, rv, &shift_minus, type,  dF,f);
			}
		}
	}
	for (i=0;i<nion;i++) 
	{
		int j; 
		j = sort_index[i]; 
		fx[j] += f[i].x;
		fy[j] += f[i].y;
		fz[j] += f[i].z;
	}
	energyInfo->sion.xx += sion_xx; 
	energyInfo->sion.xy += sion_xy; 
	energyInfo->sion.xz += sion_xz; 
	energyInfo->sion.yy += sion_yy; 
	energyInfo->sion.yz += sion_yz; 
	energyInfo->sion.zz += sion_zz; 
	energyInfo->eion+= pair_energy; 
//	printf("%d: e=%f\n",getRank(0),pair_energy); 
	heapFree(block_number); 
	HPM_Stop("pairEnergy");
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
