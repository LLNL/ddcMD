#include <string.h>
#include <limits.h>
#include "system.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "energyInfo.h"
#include "expandbuffer.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "bioGid.h"
#include "gid.h"
#include "units.h"
#include "mpiUtils.h"

int compare(const void * a, const void * b);

char *strcat_(char *a, char *b)
{
    return strcat(a, b);
}


void connectiveEnergy(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms, ETYPE *e);

int resNameDiffer(char* resName, char *name)
{
    if (strcmp(resName, name) == 0) return 0; //resClean
    return 1;
}

int countResidues(int n, GID_ORDER *gidOrder)
{
    int listSize = 0;
    gid_type gidLast = invalidGid;
    for (int j = 0; j < n; j++) //count residues
    {
        gid_type gid = gidOrder[j].gid;
        if ((gid & molResMask) != (gidLast & molResMask))
        {
            listSize++;
            gidLast = gid;
        }
    }
    return listSize;
}

void charmmResidues(SYSTEM*sys, CHARMMPOT_PARMS *parms)
{
    int incr = 4096;
    //unsigned nion=sys->nion; 
    unsigned nn = sys->nlocal;
    STATE* state = sys->collection->state;
    SETLIST *residueSet = &parms->residueSet;
    parms->gidOrder = ExpandBuffers((void*) parms->gidOrder, sizeof (GID_ORDER), nn, incr, LOCATION("charmmResidues"), "parms->gidOrder");
    GID_ORDER *gidOrder = parms->gidOrder;
    for (unsigned i = 0; i < nn; i++)
    {
        gidOrder[i].id = i;
        gidOrder[i].gid = state->label[i];
    }
    qsort(gidOrder, nn, sizeof (GID_ORDER), compare);

    int nResidue = countResidues(nn, gidOrder);
    residueSet->list = ExpandBuffers((void*) residueSet->list, sizeof (LISTNODE), nResidue, incr, LOCATION("charmmConvalent"), "residueSet->list");
    LISTNODE* residueList = residueSet->list;
    int molSize = 0;
    int listSize = 0;
    gid_type gidLast = invalidGid;
    for (unsigned j = 0; j < nn; j++) //count residues
    {
        unsigned i = gidOrder[j].id;
        gid_type gid = gidOrder[j].gid;
        if ((gid & molResMask) != (gidLast & molResMask))
        {
            RESI_CONN * resiConn = findResiConnNew(parms->charmmParms, state->species[i]->name);
            int residueAtmSize = resiConn->atomListSize;
            residueList[listSize].resiStart = (gid & molResMask);
            residueList[listSize].size = residueAtmSize;
            residueList[listSize].name = state->species[i]->name;
            residueList[listSize].resiConn = resiConn;
            residueList[listSize].resID = (gid & molResMask) >> 16;
            residueList[listSize].first = j;
            if (listSize > 0) residueList[listSize - 1].cnt = j - residueList[listSize - 1].first;
            molSize += residueAtmSize;
            listSize++;
            gidLast = gid;
        }
    }
    residueList[listSize - 1].cnt = nn - residueList[listSize - 1].first;
    residueSet->molSize = molSize;
    residueSet->listSize = listSize;
}

void charmmConvalent(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    profile(CHARMM_COVALENT, START);
    charmmResidues(sys, parms);
    SETLIST *residueSet = &parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    int incr = 4096;
    STATE* state = sys->collection->state;
    unsigned nion = sys->nion;


    GID_ORDER gidOrder2[residueSet->molSize];
    GID_ORDER *gidOrder = parms->gidOrder;

    unsigned partialNion = residueSet->molSize;

    STATE *statechpad = &(parms->statechpad);
    statechpad->nlocal = state->nlocal;
    statechpad->nion = partialNion;

    static unsigned partialNionMax = 0;
    if (partialNion > partialNionMax)
    {
        statechpad->rx = ExpandBuffers((void*) statechpad->rx, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->rx");
        statechpad->ry = ExpandBuffers((void*) statechpad->ry, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->ry");
        statechpad->rz = ExpandBuffers((void*) statechpad->rz, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->rz");

        statechpad->vx = ExpandBuffers((void*) statechpad->vx, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->vx");
        statechpad->vy = ExpandBuffers((void*) statechpad->vy, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->vy");
        statechpad->vz = ExpandBuffers((void*) statechpad->vz, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->vz");

        statechpad->fx = ExpandBuffers((void*) statechpad->fx, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->fx");
        statechpad->fy = ExpandBuffers((void*) statechpad->fy, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->fy");
        statechpad->fz = ExpandBuffers((void*) statechpad->fz, sizeof (double), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->fz");

        statechpad->label = ExpandBuffers((void*) statechpad->label, sizeof (gid_type), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->label");
        statechpad->species = ExpandBuffers((void*) statechpad->species, sizeof (SPECIES*), partialNion, incr, LOCATION("charmmConvalent"), "statechpad->species");
        partialNionMax = partialNion;
    }

    unsigned s_id = 0;
    int sp_id = 0;
    for (int i = nion; i < residueSet->molSize; i++)
    {
        gidOrder2[i].id = INT_MAX; //this will be removed in next commit. see note below on why we set gid to this high number. 
        statechpad->label[i] = INT_MAX;
        ;
    }
    unsigned index = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN * resiConn = residueList[i].resiConn;
        char * name = resiConn->resName;
        char *terminus = "x";
        if (resiConn->nTer) { terminus = "n"; }
        if (resiConn->cTer) { terminus = "c";}
        char atmName[11];
        strcpy(atmName, name);
        strcat_(atmName, terminus);
        int len = strlen(atmName);
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            //strcat_(atmName, resiConn->atomList[j]->atmName);
            strcpy(atmName + len, resiConn->atomList[j]->atmName);
            statechpad->rx[sp_id] = -1;
            s_id = gidOrder[index].id;
            if (index < sys->nion && state->species[s_id] && (state->species[s_id]->name) &&(strcmp(atmName, state->species[s_id]->name) == 0))
            {
                statechpad->label[sp_id] = state->label[s_id];
                statechpad->rx[sp_id] = state->rx[s_id];
                statechpad->ry[sp_id] = state->ry[s_id];
                statechpad->rz[sp_id] = state->rz[s_id];
                statechpad->fx[sp_id] = state->fx[s_id];
                statechpad->fy[sp_id] = state->fy[s_id];
                statechpad->fz[sp_id] = state->fz[s_id];
                statechpad->species[sp_id] = state->species[s_id];
                gidOrder2[sp_id].id = gidOrder[index].id;
                index++;
            }
            else
            { //we have a hole in the padded array, some residue atoms are missing
                //set gidOrder2 to
                //unrealistically high value: we use this info in energy functions so we know when to ignore terms
                //but we can just use the exists array mentioned below
                //this will be removed in next commit
                statechpad->label[sp_id] = INT_MAX;
                gidOrder2[sp_id].id = INT_MAX;
            }
            sp_id++;
        }
    }

    //RESRANGE resRangeTmp[count];
    RESRANGE resRange2[residueSet->listSize];
    unsigned start = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        resRange2[i].start = start;
        start += residueList[i].size;
        resRange2[i].end = start;
    }

    /*
        //read weights?
        if(parms->weightedCharmm){
          if(parms->charmmParms->charmmWeights->useFiles)
            updateWeights(parms->charmmParms->charmmWeights);
        } 
     */

    BIOENERGIES bioEnergiesZero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    parms->bioEnergies = bioEnergiesZero;

    profile(CHARMM_CONNECTIVE, START);
    for (int i = 0; i < residueSet->listSize; i++)
    {
        connectiveEnergy(statechpad, gidOrder2, sys->nlocal, residueSet->molSize, residueList[i].name, &resRange2[i], parms, e);
    }
    double etot = 0.0;
    etot += parms->bioEnergies.bond;
    etot += parms->bioEnergies.angle;
    etot += parms->bioEnergies.ub;
    etot += parms->bioEnergies.torsion;
    etot += parms->bioEnergies.impr;
    etot += parms->bioEnergies.cmap;
    etot += parms->bioEnergies.bpair;

    e->eion += etot;


    profile(CHARMM_CONNECTIVE, END);

    // Update the force to ddc items 
    // May not need to update coordinates.
    // transfer state from padded sorted state array to system state array
    index = 0;
    for (int i = 0; i < residueSet->molSize; i++)
    {
        //if (statechpad->label[i] == INT_MAX) {continue;};
        //if (gidOrder2[i].id == INT_MAX) {continue;};
        unsigned index = gidOrder2[i].id;
        if (index > sys->nion)
        {
            continue;
        }
        state->rx[index] = statechpad->rx[i];
        state->ry[index] = statechpad->ry[i];
        state->rz[index] = statechpad->rz[i];
        state->fx[index] = statechpad->fx[i];
        state->fy[index] = statechpad->fy[i];
        state->fz[index] = statechpad->fz[i];
        //index++; 
    }
    //printf("cpu bond sum %f\n", parms->bioEnergies.bond);
    profile(CHARMM_COVALENT, END);

}
