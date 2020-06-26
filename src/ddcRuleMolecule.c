#include <string.h>
#include <assert.h>

#include "simulate.h"
#include "molecule.h"
#include "system.h"
#include "gid.h"
#include "ddcRule.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "utilities.h"
#include "mpiUtils.h"
#include "heap.h"
#include "alltoall_sparse.h"

extern FILE *ddcfile;

typedef struct moleSorting_struct
{
    gid_type gid;
    unsigned int index, speciesTypeIndex, distRank, origRank;
} MOLESORTING;

void ddcRuleMolecule_init(DDCRULE *ddcRule)
{
    //DDC* ddc = (DDC*) ddcRule->parent;
    //SIMULATE* simulate = (SIMULATE*) (ddc->parent);
    //SYSTEM* system = simulate->system;
}

integer getDest(const void *element)
{
    const MOLESORTING *item = (const MOLESORTING *) element;
    return ((item->gid >> 32) % getSize(0));
}

integer getOrig(const void *element)
{
    const MOLESORTING *item = (const MOLESORTING *) element;
    return item->origRank;
}

void ddcRuleMolecule(DDC* ddc)
{
    //DDCRULE* ddcRule = ddc->ddcRule;
    SYSTEM *sys = system_getSystem(NULL);
    //SPECIES **species = sys->species; 
    PARTICLE *particles = ddc->particles;
    int *speciesIndexToMoleculeIndex = sys->moleculeClass->speciesIndexToMoleculeIndex;
    MOLECULETYPE **moleculeTypes = sys->moleculeClass->moleculeTypes;
    unsigned origRank = getRank(0);

    unsigned listBlk;
    MOLESORTING *list = (MOLESORTING *) heapGet(&listBlk);
    int listSizeMax = heapAvailable() / sizeof (*list);
    int nSend = 0;
    for (unsigned i = 0; i < ddc->number_local; i++)
    {
        int speciesTypeIndex = (particles[i].type & 0xffff);
        int moleTypeIndex = speciesIndexToMoleculeIndex[speciesTypeIndex];
        //MOLECULETYPE *moleuleType=moleculeTypes[moleTypeIndex]; 
        gid_type gid = particles[i].global_index;
        if (moleculeTypes[moleTypeIndex]->nSpecies > 1)
        {

            list[nSend].gid = gid;
            list[nSend].index = i;
            list[nSend].speciesTypeIndex = speciesTypeIndex;
            list[nSend].distRank = particles[i].domain_id;
            list[nSend].origRank = origRank;
            //      fprintf(ddcfile,"Send %16.16llx %d %d %d %d\n",list[nSend].gid,list[nSend].index,list[nSend].speciesTypeIndex,list[nSend].distRank,list[nSend].origRank); 
            nSend++;
        }
    }
    heapEndBlock(listBlk, 4 * nSend * sizeof (*list));
    //integer alltoall_sparse(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm, integer nloc,integer nlocmax, integer sz,void *data,integer (*get_dest)(const void *element));
    int nRecv = alltoall_sparse(0, getSize(0), origRank, 1, COMM_LOCAL, nSend, listSizeMax, sizeof (*list), (void*) list, getDest);
    //  for (int i=0;i<nRecv;i++) fprintf(ddcfile,"Recv %16.16llx %d %d %d %d\n",list[i].gid,list[i].index,list[i].speciesTypeIndex,list[i].distRank,list[i].origRank); 
    qsort(list, nRecv, sizeof (*list), (int(*)(const void*, const void*))compareGid);
    //   for (int i=0;i<nRecv;i++) fprintf(ddcfile,"Sort %16.16llx %d %d %d %d\n",list[i].gid,list[i].index,list[i].speciesTypeIndex,list[i].distRank,list[i].origRank); 
    int cnt = 0;
    for (int i = 0; i < nRecv; i++)
    {
        int distRank = 0; //TODO: distRank has to given a initial value
        //int nSpecies;
        if ((list[i].gid & 0xffffffff) == 0)
        {
            assert(cnt == 0);
            int speciesTypeIndex = list[i].speciesTypeIndex;
            int moleTypeIndex = speciesIndexToMoleculeIndex[speciesTypeIndex];
            //MOLECULETYPE *moleuleType=moleculeTypes[moleTypeIndex]; 
            int offset = moleculeTypes[moleTypeIndex]->ownershipSpeciesOffset;
            cnt = moleculeTypes[moleTypeIndex]->nSpecies;
            distRank = list[i + offset].distRank;
        }
        list[i].distRank = distRank;
        cnt--;
        //fprintf(ddcfile,"Sorted %16.16llx %d %d %d %d\n",list[i].gid,list[i].index,list[i].speciesTypeIndex,list[i].distRank,list[i].origRank); 
    }
    assert(cnt == 0);
    nSend = nRecv;
    nRecv = alltoall_sparse(0, getSize(0), origRank, 1, COMM_LOCAL, nSend, listSizeMax, sizeof (*list), (void*) list, getOrig);
    heapFree(listBlk);
    for (int i = 0; i < nRecv; i++) particles[list[i].index].domain_id = list[i].distRank;
}
