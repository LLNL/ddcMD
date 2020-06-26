#include <string.h>

#include "simulate.h"
#include "system.h"
#include "gid.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioGid.h"
#include "ddcRule.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "utilities.h"
#include "mpiUtils.h"
#include <assert.h>

typedef struct ddcResidue_st
{
    int start, size;
    RESI_CONN* residue;
} DDCRESIDUE;
extern FILE *ddcfile;

void charmmResidues(SYSTEM *sys, CHARMMPOT_PARMS *parms);
int sortParticlesByGid(const PARTICLE*p1, const PARTICLE*p2);

void ddcRuleMartini_init(DDCRULE *ddcRule)
{
    DDC* ddc = (DDC*) ddcRule->parent;
    SIMULATE* simulate = (SIMULATE*) (ddc->parent);
    SYSTEM* system = simulate->system;
    int foundMMFF = 0;
    for (int i = 0; i < system->npotential; ++i)
    {
        if (system->potential[i]->itype == MARTINI)
        {
            CHARMMPOT_PARMS* charmmpot_parms = (CHARMMPOT_PARMS *) (system->potential[i]->parms);
            ddcRule->parms = charmmpot_parms;
            foundMMFF = 1;
            break;
        }
    }
    if (foundMMFF == 0)
    {
        error_action("ddcRule Error: Cannot find CHARMM_PARMS in the system potentials", ERROR_IN("ddcRule", ABORT));
    }
}

int completeResidueError(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    SETLIST residueSet = parms->residueSet;
    charmmResidues(sys, parms);
    int flag = 0;
    for (int i = 0; i < residueSet.listSize; i++)
    {
        LISTNODE residue = residueSet.list[i];
        if (residue.size != residue.cnt)
        {
            flag = 1;
        }
    }
    return flag;
}

gid_type checkResidue(int size, SYSTEM *sys, PARTICLE *particles, CHARMMPOT_PARMS *parms, int *completeResidue)
{
    int speciesIndex = (particles[0].type & 0xffff);
    char *name = sys->species[speciesIndex]->name;
    RESI_CONN *resiConn = findResiConnNew(parms->charmmParms, name);
    gid_type centerAtom = resiConn->centerAtom;
    int nAtom = resiConn->atomListSize;
    *completeResidue = 0;
    if (size == nAtom) *completeResidue = 1;
    gid_type residueId = invalidGid;
    for (int k = 0; k < size; k++)
    {
        gid_type gid = particles[k].global_index;
        if ((gid & atmgrpMask) == centerAtom)
        {
            unsigned id = particles[k].domain_id;
            residueId = (gid & molresMask) + id;
            break;
        }
    }
    return residueId;
}

int ddcCountResidues(int n, PARTICLE *particles, SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    gid_type hasCenter[1024]; // Need to remove the fixed length

    struct
    {
        gid_type gid;
        int start;
        int size;
    } noCenter[1024]; // Need to remove the fixed length
    int nHasCenter = 0;
    int nNoCenter = 0;

    gid_type gidLast = invalidGid;
    int start = 0;
    int cnt = 0;
    for (int j = 0; j < n; j++)
    {
        gid_type gid = particles[j].global_index;
        if ((gid & molResMask) != (gidLast & molResMask))
        {
            int size = j - start;
            if (size > 0)
            {
                int completeResidue = 0;
                gid_type residueId = checkResidue(size, sys, particles + start, parms, &completeResidue);
                if (!completeResidue)
                {
                    if (residueId != invalidGid)
                    {
                        hasCenter[nHasCenter++] = residueId;
                    }
                    else
                    {
                        noCenter[nNoCenter].gid = (gidLast & molresMask);
                        noCenter[nNoCenter].start = start;
                        noCenter[nNoCenter].size = size;
                        nNoCenter++;
                    }
                }
            }
            start = j;
            gidLast = gid;
            cnt++;
        }
    }
    int size = n - start;
    if (size > 0)
    {
        int completeResidue = 0;
        gid_type residueId = checkResidue(size, sys, particles + start, parms, &completeResidue);
        if (!completeResidue)
        {
            if (residueId != invalidGid)
            {
                hasCenter[nHasCenter++] = residueId;
                for (int i = 0; i < size; i++) particles[i + start].domain_id = (residueId & 0xffff);
            }
            else
            {
                noCenter[nNoCenter].gid = (gidLast & molresMask);
                noCenter[nNoCenter].start = start;
                noCenter[nNoCenter].size = size;
                nNoCenter++;
            }
        }
    }
    int nSize = getSize(0);
    int ncArray[nSize], disp[nSize + 1];

    MPI_Allgather(&nHasCenter, 1, MPI_INT, ncArray, 1, MPI_INT, COMM_LOCAL); //Gather nHasCenter; 

    disp[0] = 0;
    for (int i = 1; i <= nSize; i++) disp[i] = disp[i - 1] + ncArray[i - 1]; //find displacments. 
    if (disp[nSize] == 0) return cnt;
    gid_type centers[disp[nSize]];

    MPI_Allgatherv(hasCenter, nHasCenter, MPI_GID_TYPE, centers, ncArray, disp, MPI_GID_TYPE, COMM_LOCAL);

    qsort(centers, disp[nSize], sizeof (gid_type), (int(*)(const void*, const void*))compareGid);
    qsort(noCenter, nNoCenter, sizeof (noCenter[0]), (int(*)(const void*, const void*))compareGid);
    int index = 0;
    for (int i = 0; i < nNoCenter; i++)
    {
        for (; index < disp[nSize]; index++)
            if ((centers[index] & molresMask) == (noCenter[i].gid & molresMask)) break;
        assert(index != disp[nSize]);
        //      fprintf(ddcfile,"Center %4d %4d %016"PRIx64" %016"PRIx64" %4d %4d\n",i,index,centers[index],noCenter[i].gid,noCenter[i].start,noCenter[i].size); 
        for (int j = 0; j < noCenter[i].size; j++)
        {
            particles[noCenter[i].start + j].domain_id = (centers[index]&0xffff);
            //          fprintf(ddcfile,"Center %4d %4d %4d %016"PRIx64" %016"PRIx64"\n",i,j,index,centers[index],particles[noCenter[i].start+j].global_index);
        }
    }
    //   for (int i=0;i<disp[nSize];i++) fprintf(ddcfile,"centers %d %016"PRIx64"\n",i,centers[i]); 
    return cnt;
}
// misc debugging code 
//fprintf(ddcfile,"x %4d %4d %016"PRIx64" %016"PRIx64" %4d %4d\n",i,index,centers[index],noCenter[i].gid,noCenter[i].start,noCenter[i].size); 
//for (int i=0;i<nSize;i++) fprintf(ddcfile,"nc disp %d %d %d\n",i,ncArray[i],disp[i]); 
//fprintf(ddcfile,"%d: nNoCenter=%d nHasCenter=%d nAll=%d\n",getRank(0),nNoCenter,nHasCenter,disp[nSize]); 
//for (int i=0;i<nHasCenter;i++) fprintf(ddcfile,"hasCenter %d %016"PRIx64"\n",i,hasCenter[i]); 
//for (int i=0;i<disp[nSize];i++) fprintf(ddcfile,"BEFORE SORT %d %016"PRIx64"\n",i,centers[i]); 
//for (int i=0;i<disp[nSize];i++) fprintf(ddcfile,"AFTER SORT %d %016"PRIx64"\n",i,centers[i]); 
//for (int i=0;i<nNoCenter;i++) fprintf(ddcfile,"SORTED noCenter %d %016"PRIx64"\n",i,noCenter[i].gid); 

void ddcRuleMartini(DDC* ddc)
{

    DDCRULE* ddcRule = ddc->ddcRule;
    CHARMMPOT_PARMS *parms = (CHARMMPOT_PARMS *) ddcRule->parms;
    SYSTEM *sys = system_getSystem(NULL);
    for (unsigned i = 0; i < ddc->number_local; i++) ddc->particles[i].ifirst = i; /*Use ifirst as temp storage for the sort order*/
    qsort(ddc->particles, ddc->number_local, sizeof (PARTICLE), (int(*)(const void*, const void*))sortParticlesByGid); /*Sorts ddc copy of particles info by globalindex */
    //   for (int ii=0;ii<ddc->number_local;ii++)
    //  fprintf(ddcfile,"Before %5d %016"PRIx64" %2d\n"   , ii,ddc->particles[ii].global_index,ddc->particles[ii].domain_id); 
    ddcCountResidues(ddc->number_local, ddc->particles, sys, parms);
    //   for (int ii=0;ii<ddc->number_local;ii++)
    //  fprintf(ddcfile,"After %5d %016"PRIx64" %2d\n"   , ii,ddc->particles[ii].global_index,ddc->particles[ii].domain_id); 
}
