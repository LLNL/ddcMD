#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "three_algebra.h"
#include "object.h"
#include "pio.h"
/*#include "error.h"*/
#include "utilities.h"
#include "crc32.h"
#include "simulate.h"
#include "system.h"
#include "ddcMalloc.h"
#include "io.h"
#include "box.h"
#include "units.h"
#include "uthash.h"
#include "ioUtils.h"
#include "preduce.h"

enum SMEAR_METHOD
{
	SMEAR_IMPULSE, SMEAR_HAT
};


/** For proper hashing the label and species members *must* be
 * adjacent in this struct. */
typedef struct grid_point_st
{
   LONG64 label;
   unsigned species;
   double nParticles;
   THREE_VECTOR p,K,E;
   double mass, U, virial, ESpotential;
	THREE_SMATRIX virialTensor;
   UT_hash_handle hh; // make structure hash-able
} GRID_POINT;

typedef struct coarsegrain_parms_st
{
	int outputMode;
   int nx, ny, nz;
   char *filename; 
   int nfiles;
   char *misc_info;
   char **field_units;
   int nfields; 
	int nsteps;
	double smearRadius;
	enum SMEAR_METHOD smearMethod;
	unsigned keylen;
	GRID_POINT* hashHead;
}  COARSEGRAIN_PARMS ; 

typedef struct cg_helper_st
{
	unsigned ny, nz;
	unsigned keylen;
	THREE_VECTOR deltai;
	THREE_VECTOR scaled_corner;
	LONG64* lastLabel;
	GRID_POINT** lastCell;
	GRID_POINT** hashHead;
} CG_HELPER;

static THREE_SMATRIX sm_zero = {0.0, 0.0, 0.0,   0.0, 0.0, 0.0};

int getRank(int); 

static void cg_addPairGuts(
	CG_HELPER* helper, LONG64 label, int iSpecies,
	THREE_VECTOR* r12, THREE_VECTOR* f12, double v12, double weight);

COARSEGRAIN_PARMS *coarsegrain_parms(ANALYSIS *analysis)
{
	char string[1024]; 
	char* smearString; 
	SIMULATE* simulate = (SIMULATE *)analysis->parent; 
	SYSTEM* sys=simulate->system; 
	SPECIES** sp = sys->species; 
	int nspecies = sys->nspecies; 

	COARSEGRAIN_PARMS* parms = analysis->parms= (COARSEGRAIN_PARMS *) ddcMalloc(sizeof(COARSEGRAIN_PARMS));
	OBJECT* obj = (OBJECT*) analysis;
	object_get(obj, "nx", &parms->nx, INT, 1, "0");
	object_get(obj, "ny", &parms->ny, INT, 1, "0");
	object_get(obj, "nz", &parms->nz, INT, 1, "0");
	object_get(obj, "filename", &parms->filename, STRING, 1, "cgrid");
	object_get(obj, "outputMode", &parms->outputMode, INT, 1, "2");
	object_get(obj, "nfiles", &parms->nfiles, INT, 1, "0");
	object_get(obj, "smearRadius", &parms->smearRadius, WITH_UNITS, 1, "0", "l", NULL);
	parms->nfields = object_getv(obj, "field_units", (void*)&parms->field_units, STRING,ABORT_IF_NOT_FOUND);
	object_get(obj, "smearMethod", &smearString, STRING, 1, "impulse");
	parms->smearMethod = SMEAR_IMPULSE;
	if (strcasecmp(smearString, "hat") == 0)
		parms->smearMethod = SMEAR_HAT;
	parms->nsteps = 0;
	parms->keylen = offsetof(GRID_POINT, species) - offsetof(GRID_POINT, label) + sizeof(parms->hashHead->species);
	parms->hashHead = NULL;
	
	sprintf(string,"nx = %d; ny = %d; nz=%d; smearRadius=%f; smearMethod=%s;\n"
			  "       ouputrate=%d; eval_rate=%d; species=",
			  parms->nx,parms->ny,parms->nz, units_convert(parms->smearRadius, NULL, "l"), smearString,
			  analysis->outputrate, analysis->eval_rate); 
	for (int j=0;j<nspecies;j++)
	{
		strcat(string," ");
		strcat(string,sp[j]->name);
	}
	strcat(string,";");
	parms->misc_info = strdup(string); 
	return parms; 
}

/** Allocate and initialize the CG_HELPER object that stores precomputed
 * values that can be re-used for the duration of a time step.  Helper
 * object also contains pointers to certain quantities in the
 * COARSEGRAIN_PARMS that are needed to add data to the hash tables.
 * Caller must call cg_helperFree when the helper object is no longer
 * needed.  Returns NULL if this is a time step that does not require
 * CoarseGrain evaluation. */
CG_HELPER* cg_helperInit(void)
{
	SIMULATE* simulate=simulate_getSimulate(NULL);

	ANALYSIS* cgObj = NULL;
	for (int ii=0; ii< simulate->nanalysis; ++ii)
		if (simulate->analysis[ii]->itype == COARSEGRAIN)
		{
			cgObj = simulate->analysis[ii];
			break;
		}
	if (!cgObj )
		return NULL; 

	if (! TEST0(simulate->loop, cgObj->eval_rate) )
		return NULL;

	COARSEGRAIN_PARMS* parms = (COARSEGRAIN_PARMS*) cgObj->parms;
	CG_HELPER* helper = ddcMalloc(sizeof(CG_HELPER));

	helper->hashHead = &parms->hashHead;
	helper->ny = parms->ny;
	helper->nz = parms->nz;
	helper->keylen = parms->keylen;
	
	THREE_VECTOR corner = box_get_corner(simulate->system->box);
	THREE_VECTOR bbox = box_get_boundingbox(simulate->system->box);
	helper->deltai.x = parms->nx/bbox.x; 
	helper->deltai.y = parms->ny/bbox.y; 
	helper->deltai.z = parms->nz/bbox.z; 
	helper->scaled_corner.x = corner.x * helper->deltai.x; 
	helper->scaled_corner.y = corner.y * helper->deltai.y; 
	helper->scaled_corner.z = corner.z * helper->deltai.z; 

	int nSpecies = simulate->system->nspecies; 
	helper->lastLabel = ddcMalloc(nSpecies * sizeof(LONG64));
	helper->lastCell  = ddcMalloc(nSpecies * sizeof(GRID_POINT*));
	
	for (int ii=0; ii<nSpecies; ++ii)
	{
		helper->lastLabel[ii] = parms->nx * parms->ny * parms->nz;
		helper->lastCell[ii] = NULL;
	}

	return helper;
}

void cg_helperFree(CG_HELPER* h)
{
	if (!h)
		return;
	ddcFree(h->lastCell);
	ddcFree(h->lastLabel);
	ddcFree(h);
}


void cg_addPair(CG_HELPER* helper, THREE_VECTOR* r1, int s1, THREE_VECTOR* r2, int s2, THREE_VECTOR* r12, THREE_VECTOR* f12, double v12)
{
	THREE_INT ig1,ig2;
	ig1.x=r1->x*helper->deltai.x - helper->scaled_corner.x; 
	ig1.y=r1->y*helper->deltai.y - helper->scaled_corner.y; 
	ig1.z=r1->z*helper->deltai.z - helper->scaled_corner.z; 
	LONG64 label1 = ig1.z+helper->nz*(ig1.y + helper->ny*ig1.x);
	ig2.x=r2->x*helper->deltai.x - helper->scaled_corner.x; 
	ig2.y=r2->y*helper->deltai.y - helper->scaled_corner.y; 
	ig2.z=r2->z*helper->deltai.z - helper->scaled_corner.z; 
	LONG64 label2 = ig2.z+helper->nz*(ig2.y + helper->ny*ig2.x);
	cg_addPairGuts(helper, label1, s1, r12, f12, v12, 0.5);
	cg_addPairGuts(helper, label2, s2, r12, f12, v12, 0.5);
}

static void cg_addPairGuts(
	CG_HELPER* helper, LONG64 label, int iSpecies,
	THREE_VECTOR* r12, THREE_VECTOR* f12, double v12, double weight)
{
	GRID_POINT* here=NULL;
	if (label == helper->lastLabel[iSpecies]) here = helper->lastCell[iSpecies];
	else
	{
		GRID_POINT probe;
		probe.label = label;
		probe.species = iSpecies;
		HASH_FIND(hh, *helper->hashHead, &probe.label, helper->keylen, here);
	}
	if (!here)
	{
		here = calloc(1, sizeof(GRID_POINT));
		here->label = label;
		here->species = iSpecies;
		here->nParticles = 0.0;
		here->p = vzero;
		here->K = vzero;
		here->E = vzero;
		here->mass = 0.0;
		here->U = 0.0;
		here->ESpotential = 0.0;
		here->virial = 0.0;
		here->virialTensor = sm_zero;
		HASH_ADD(hh, *helper->hashHead, label, helper->keylen, here);
	}
	helper->lastLabel[iSpecies] = label;
	helper->lastCell[iSpecies] = here;
// Note Trace(virialTensor)/3 = virial; 
	here->virialTensor.xx += weight*f12->x*r12->x;
	here->virialTensor.yy += weight*f12->y*r12->y;
	here->virialTensor.zz += weight*f12->z*r12->z;
	here->virialTensor.xy += weight*f12->x*r12->y;
	here->virialTensor.xz += weight*f12->x*r12->z;
	here->virialTensor.yz += weight*f12->y*r12->z;
	here->U += weight*v12; 
}


void coarsegrain_eval(ANALYSIS* analysis)
{
	COARSEGRAIN_PARMS* parms = (COARSEGRAIN_PARMS*) analysis->parms; 
	SIMULATE* simulate = (SIMULATE*) analysis->parent; 
	SYSTEM* sys = simulate->system; 
	STATE* state = sys->collection->state; 
	SPECIES** species = state->species; 
	unsigned nlocal = sys->nlocal; 

	double* rx = state->rx; 
	double* ry = state->ry; 
	double* rz = state->rz; 
	double* vx = state->vx; 
	double* vy = state->vy; 
	double* vz = state->vz;
	double* fx = state->fx;
	double* fy = state->fy;
	double* fz = state->fz;
	double* potential = state->potentialEnergy; 
	double* virial = state->virial; 
	THREE_SMATRIX* sion = state->sion;
	
	THREE_VECTOR corner = box_get_corner(sys->box);
	THREE_VECTOR bbox = box_get_boundingbox(NULL);
	THREE_VECTOR deltai,scaled_corner;
	deltai.x = parms->nx/bbox.x; 
	deltai.y = parms->ny/bbox.y; 
	deltai.z = parms->nz/bbox.z; 
	scaled_corner.x = corner.x * deltai.x; 
	scaled_corner.y = corner.y * deltai.y; 
	scaled_corner.z = corner.z * deltai.z; 


	THREE_VECTOR lSmearInv = vzero;
	THREE_VECTOR lSmearHalf = vzero;
	if (parms->smearRadius > 0)
	{
		THREE_VECTOR lSmear;
		lSmear.x = MIN( 2.0*parms->smearRadius, bbox.x/(1.0*parms->nx));
		lSmear.y = MIN( 2.0*parms->smearRadius, bbox.y/(1.0*parms->ny));
		lSmear.z = MIN( 2.0*parms->smearRadius, bbox.z/(1.0*parms->nz));

		lSmearInv.x = 1.0/lSmear.x;
		lSmearInv.y = 1.0/lSmear.y;
		lSmearInv.z = 1.0/lSmear.z;

		lSmearHalf.x = 0.5 * lSmear.x;
		lSmearHalf.y = 0.5 * lSmear.y;
		lSmearHalf.z = 0.5 * lSmear.z;
	}
	
	parms->nsteps++;
	for (unsigned iAtom=0; iAtom<nlocal; ++iAtom) 
	{
		double rxi = rx[iAtom];
		double ryi = ry[iAtom];
		double rzi = rz[iAtom];
		backInBox_fast(&rxi, &ryi, &rzi);
		THREE_VECTOR r;
		r.x=rxi*deltai.x - scaled_corner.x; 
		r.y=ryi*deltai.y - scaled_corner.y; 
		r.z=rzi*deltai.z - scaled_corner.z; 

		int spread = 2;
		int igx[2], igy[2], igz[2];
		double wx[2], wy[2], wz[2];
		if ( parms->smearRadius <= 0)
		{
			spread = 1;
			igx[0] = r.x; igy[0] = r.y; igz[0] = r.z;
			wx[0]  = 1.0; wy[0]  = 1.0; wz[0]  = 1.0;
		}
		else
		{
			spread = 2;
			THREE_INT iWall;
			iWall.x = floor(r.x + 0.5);
			iWall.y = floor(r.y + 0.5);
			iWall.z = floor(r.z + 0.5);
			
			THREE_VECTOR delta;
			delta.x = iWall.x - r.x;
			delta.y = iWall.y - r.y;
			delta.z = iWall.z - r.z;

			delta.x = MIN(delta.x, lSmearHalf.x);
			delta.x = MAX(delta.x, -lSmearHalf.x);
			delta.y = MIN(delta.y, lSmearHalf.y);
			delta.y = MAX(delta.y, -lSmearHalf.y);
			delta.z = MIN(delta.z, lSmearHalf.z);
			delta.z = MAX(delta.z, -lSmearHalf.z);
			
			igx[0] = iWall.x-1; if (igx[0] == -1)        igx[0] = parms->nx-1; 
			igy[0] = iWall.y-1; if (igy[0] == -1)        igy[0] = parms->ny-1; 
			igz[0] = iWall.z-1; if (igz[0] == -1)        igz[0] = parms->nz-1;
			igx[1] = iWall.x;   if (igx[1] == parms->nx) igx[1] = 0;
			igy[1] = iWall.y;   if (igy[1] == parms->ny) igy[1] = 0;
			igz[1] = iWall.z;   if (igz[1] == parms->nz) igz[1] = 0;
			
			switch (parms->smearMethod)
			{
			  case SMEAR_IMPULSE:
				wx[0] = 0.5 + (delta.x * lSmearInv.x);
				wy[0] = 0.5 + (delta.y * lSmearInv.y);
				wz[0] = 0.5 + (delta.z * lSmearInv.z);
				break;
			  case SMEAR_HAT:
				wx[0] = 0.5 + 2*delta.x*lSmearInv.x*(1.0 - fabs(delta.x)*lSmearInv.x);
				wy[0] = 0.5 + 2*delta.y*lSmearInv.y*(1.0 - fabs(delta.y)*lSmearInv.y);
				wz[0] = 0.5 + 2*delta.z*lSmearInv.z*(1.0 - fabs(delta.z)*lSmearInv.z);
				break;
			  default:
				assert (1==0);
			}
			

			assert(wx[0] >= 0);
			assert(wy[0] >= 0);
			assert(wz[0] >= 0);
			assert(wx[0] <= 1);
			assert(wy[0] <= 1);
			assert(wz[0] <= 1);
			wx[1] = 1.0 - wx[0];
			wy[1] = 1.0 - wy[0];
			wz[1] = 1.0 - wz[0];
		}
		
		const unsigned iSpecies = species[iAtom]->index; 
		const double mi = ((ATOMTYPE_PARMS*)(species[iAtom])->parm)->mass;
		const double qi = ((ATOMTYPE_PARMS*)(species[iAtom])->parm)->charge;
		const double px = mi*vx[iAtom];
		const double py = mi*vy[iAtom];
		const double pz = mi*vz[iAtom];
		const double Kx = 0.5*mi*(vx[iAtom]*vx[iAtom]);
		const double Ky = 0.5*mi*(vy[iAtom]*vy[iAtom]);
		const double Kz = 0.5*mi*(vz[iAtom]*vz[iAtom]);
		const double Ex = fx[iAtom]/qi;
		const double Ey = fy[iAtom]/qi;
		const double Ez = fz[iAtom]/qi;
		const double U = potential[iAtom];
		const double ESpotential = U/qi;
		const double vir = virial[iAtom];
		const double vir_xx = sion[iAtom].xx;
		const double vir_yy = sion[iAtom].yy;
		const double vir_zz = sion[iAtom].zz;
		const double vir_xy = sion[iAtom].xy;
		const double vir_xz = sion[iAtom].xz;
		const double vir_yz = sion[iAtom].yz;
		for (int ii=0; ii<spread; ++ii)
		{
			for (int jj=0; jj<spread; ++jj)
			{
				for (int kk=0; kk<spread; ++kk)
				{
					double weight = wx[ii]*wy[jj]*wz[kk];
					if (weight < 1e-20)
						continue;
					
					
					unsigned label = igz[kk]+parms->nz*(igy[jj] + parms->ny*igx[ii]);
					GRID_POINT probe;
					probe.label = label;
					probe.species = iSpecies;
					GRID_POINT* here = NULL;
					HASH_FIND(hh, parms->hashHead, &probe.label, parms->keylen, here);
					if (!here)
					{
						here = calloc(1, sizeof(GRID_POINT));
						here->label = label;
						here->species = iSpecies;
						here->nParticles = 0;
						here->p = vzero;
						here->K = vzero;
						here->E = vzero;
						here->mass = 0.0;
						here->U = 0.0;
						here->virial = 0.0;
						here->ESpotential = 0.0;
						HASH_ADD(hh, parms->hashHead, label, parms->keylen, here);
					}
					
					here->p.x        += weight * px;
					here->p.y        += weight * py;
					here->p.z        += weight * pz;
					here->K.x        += weight * Kx;
					here->K.y        += weight * Ky;
					here->K.z        += weight * Kz;
					here->E.x        += weight * Ex;
					here->E.y        += weight * Ey;
					here->E.z        += weight * Ez;
					here->U          += weight * U;
					here->ESpotential+= weight * ESpotential;
					here->virial     += weight * vir;
					here->mass       += weight * mi; 
					here->nParticles += weight; 
					here->virialTensor.xx -= weight * vir_xx;
					here->virialTensor.yy -= weight * vir_yy;
					here->virialTensor.zz -= weight * vir_zz;
					here->virialTensor.xy -= weight * vir_xy;
					here->virialTensor.xz -= weight * vir_xz;
					here->virialTensor.yz -= weight * vir_yz;
				}
			}
		}
		
	}
}
void coarsegrain_output(ANALYSIS *analysis)
{
	char filename[1024];
	COARSEGRAIN_PARMS* parms = (COARSEGRAIN_PARMS *)analysis->parms; 
	SIMULATE* simulate =(SIMULATE *)(analysis->parent); 
	SYSTEM* sys=simulate->system; 
	CreateSnapshotdir(simulate, NULL);
	snprintf(filename, 1024,"%s/%s", simulate->snapshotdir,parms->filename);
	PFILE* file = Popen(filename, "w", COMM_LOCAL);
	int nspecies = sys->nspecies; 
	int nFields = 13;

	unsigned long long maxLabel = (parms->nx * parms->ny * parms->nz) - 1;
	unsigned long long maxSpecies = nspecies - 1;
	unsigned labelFieldSize = bFieldSize(maxLabel);
	unsigned speciesFieldSize = bFieldSize(maxSpecies);
	
	char fieldNames[1024], fieldTypes[1024], fieldUnits[1024];
	int lrec = 0;
	int mode = parms->outputMode;
	switch (mode)
	{
	  case 1:
		nFields = 13;
		lrec = 11*4 + labelFieldSize + speciesFieldSize;
		sprintf(fieldNames, "checksum label species_index number_particles mass  Kx Ky Kz U virial px py pz");
		sprintf(fieldTypes, "u4 b%1u b%1u f4 f4 f4 f4 f4 f4 f4 f4 f4 f4", labelFieldSize, speciesFieldSize);
		sprintf(fieldUnits, "1 1 1 1 amu eV eV eV eV eV amu*Angstrom/fs amu*Angstrom/fs amu*Angstrom/fs");
		break;
	  case 2:
		nFields = 19;
		lrec = 17*4 + labelFieldSize + speciesFieldSize;
		sprintf(fieldNames, "checksum label species_index number_particles mass  Kx Ky Kz U virial px py pz vir_xx vir_yy vir_zz vir_xy vir_xz vir_yz");
		sprintf(fieldTypes, "u4 b%1u b%1u f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4 f4", labelFieldSize, speciesFieldSize);
		sprintf(fieldUnits, "1 1 1 1 amu eV eV eV eV eV amu*Angstrom/fs amu*Angstrom/fs amu*Angstrom/fs eV eV eV eV eV eV");
		break;
	  case 3:
		nFields = 12;
		lrec = 10*4 + labelFieldSize + speciesFieldSize;
		sprintf(fieldNames, "checksum label species_index number_particles mass  px py pz Ex Ey Ez ESpotential");
		sprintf(fieldTypes, "u4 b%1u b%1u f4 f4 f4 f4 f4 f4 f4 f4 f4", labelFieldSize, speciesFieldSize);
		sprintf(fieldUnits, "1 1 1 1 amu amu*Angstrom/fs amu*Angstrom/fs amu*Angstrom/fs eV/e/Angstrom eV/e/Angstrom eV/e/Angstrom eV/e");
		break;
		
	  default:
		assert(0==1);
	}

	PioSet(file, "recordLength", lrec);
	PioSet(file, "datatype", FIXRECORDBINARY);
	PioSet(file, "numberRecords", 0llu);
	PioSet(file, "checksum", CRC32);
	PioSet(file, "nfields", nFields);
	PioSet(file, "field_names", fieldNames);
	PioSet(file, "field_types", fieldTypes);
	PioSet(file, "field_units", fieldUnits);

	unsigned tmpLen = strlen(parms->misc_info) + 20;
	char string[tmpLen];
	sprintf(string, "%s\n  nSamples = %d;", parms->misc_info, parms->nsteps);
	PioSet(file, "misc_info", string);

	PioReserve(file, lrec*parms->nx*parms->ny*parms->nz*nspecies/file->size);
	if (parms->nfiles > 0)
	   PioSet(file, "ngroup", parms->nfiles);
	if (getRank(0) == 0)
		write_fileheader(file, simulate, "cgrid" );

	double mass_convert = units_convert(1.0,NULL,"amu");
	double energy_convert = units_convert(1.0,NULL,"eV");
	double momentum_convert = units_convert(1.0,NULL,"amu*Angstrom/fs");
	double E_convert = units_convert(1.0,NULL,"eV/e/Angstrom");
	double ESpotential_convert = units_convert(1.0,NULL,"eV/e");
	double nstepsInv = 1.0/(parms->nsteps*1.0);
	unsigned char* line = (unsigned char*) ddcMalloc(lrec*sizeof(char));
	unsigned* offset = ddcMalloc((nFields+1)*sizeof(unsigned));
	parseFieldSizes(offset, fieldTypes);
	
	for (GRID_POINT* iter=parms->hashHead; iter!=NULL; iter=iter->hh.next)
   {
		//TODO: Generalize this...it's getting rather messy
            
		unsigned i4;
		LONG64   i8;
		float    f4;
		unsigned ii = 3;
		i8 = iter->label;   bFieldPack(line+offset[1], labelFieldSize,   i8);
		i8 = iter->species; bFieldPack(line+offset[2], speciesFieldSize, i8);
		f4 = iter->nParticles*nstepsInv;            copyBytes(line+offset[ii++],  &f4, 4);
		f4 = iter->mass*mass_convert*nstepsInv;     copyBytes(line+offset[ii++],  &f4, 4);
		if (mode == 1 || mode == 2)
		{
			f4 = iter->K.x*energy_convert*nstepsInv;    copyBytes(line+offset[ii++],  &f4, 4);
			f4 = iter->K.y*energy_convert*nstepsInv;    copyBytes(line+offset[ii++],  &f4, 4);
			f4 = iter->K.z*energy_convert*nstepsInv;    copyBytes(line+offset[ii++],  &f4, 4);
			f4 = iter->U*energy_convert*nstepsInv;      copyBytes(line+offset[ii++],  &f4, 4);
			f4 = iter->virial*energy_convert*nstepsInv; copyBytes(line+offset[ii++],  &f4, 4);
		}

		f4 = iter->p.x*momentum_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
		f4 = iter->p.y*momentum_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
		f4 = iter->p.z*momentum_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
		
		if (mode == 2)
		{
			f4 = iter->virialTensor.xx*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->virialTensor.yy*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->virialTensor.zz*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->virialTensor.xy*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->virialTensor.xz*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->virialTensor.yz*energy_convert*nstepsInv;  copyBytes(line+offset[ii++], &f4, 4);
		}
		else if (mode==3)
		{
			f4 = iter->E.x*E_convert*nstepsInv;                   copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->E.y*E_convert*nstepsInv;                   copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->E.z*E_convert*nstepsInv;                   copyBytes(line+offset[ii++], &f4, 4);
			f4 = iter->ESpotential*ESpotential_convert*nstepsInv; copyBytes(line+offset[ii++], &f4, 4);
		}
		
		i4 =  checksum_crc32_table(line+offset[1],lrec-offset[1]);
	   copyBytes(line, &i4, 4);
	   Pwrite(line, lrec, 1, (PFILE*) file);
	}
	Pclose(file);
	parms->nsteps = 0;
	while (parms->hashHead)
	{
		GRID_POINT* tmp = parms->hashHead;
		HASH_DEL(parms->hashHead, tmp);
		free(tmp);
	}
	ddcFree(offset);
	ddcFree(line);
}

void coarsegrain_clear(ANALYSIS *analysis)
{
	COARSEGRAIN_PARMS* parms = (COARSEGRAIN_PARMS *)analysis->parms; 
	while (parms->hashHead)
	{
		GRID_POINT* tmp = parms->hashHead;
		HASH_DEL(parms->hashHead, tmp);
		free(tmp);
	}
	parms->nsteps=0;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
