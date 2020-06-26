#include "bioMartini.h"
#include "object.h"
#include "ddcMalloc.h"
#include "utilities.h"
#include "error.h"
#include "codata.h"
#include "units.h"
#include "ptiming.h"
#include "preduce.h"
#include "external.h"
#include "mpiUtils.h"
#include "simulate.h"
#include "species.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "bioMMFF.h"
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "molecule.h"
#include "gid.h"
#include "bioGid.h"
#include "accelerator.h"
#include "HAVEGPU.h"

RESI_CONN* findResiConnNew(CHARMM_PARMS* charmmParms, char* name);
void reOrgPairs(SYSTEM*sys, CHARMMPOT_PARMS *parms);
int completeResidueError(SYSTEM *sys, CHARMMPOT_PARMS *parms);

int getMartiniLJpairParm(MMFF* mmff, char* atmTypeI, char* atmTypeJ, double* eps, double* sigma)
{
    for (int i = 0; i < mmff->nLJParms; i++)
    {
        LJPARMS *ljparms = mmff->ljParms[i];
        if (strcmp(atmTypeI, ljparms->atomTypeI) == 0 && strcmp(atmTypeJ, ljparms->atomTypeJ) == 0)
        {
            *eps = ljparms->eps;
            *sigma = ljparms->sigma;
            return 0;
        }
        if (strcmp(atmTypeI, ljparms->atomTypeJ) == 0 && strcmp(atmTypeJ, ljparms->atomTypeI) == 0)
        {
            *eps = ljparms->eps;
            *sigma = ljparms->sigma;
            return 0;
        }
    }

    return -1;
}

void validateExclusions(MMFF *mmff)
{
    for (int r = 0; r < mmff->nResiParms; ++r)
    {
        RESIPARMS* resiParm = mmff->resiParms[r];

        for (int e = 0; e < resiParm->nExclude; e++)
        {
            EXCLUDEPARMS* excludeparms = resiParm->exclusionList[e];

            for (int i = 0; i < resiParm->nBonds; ++i)
            {
                BONDPARMS* bondparms = resiParm->bondList[i];

                if (bondparms->func == 1)
                {
                    if (excludeparms->atomI == bondparms->atomI && excludeparms->atomJ == bondparms->atomJ)
                    {
                        excludeparms->valid = 0;
                        break;
                    }
                    if (excludeparms->atomI == bondparms->atomJ && excludeparms->atomJ == bondparms->atomI)
                    {
                        excludeparms->valid = 0;
                        break;
                    }
                }
            }

        }

        for (int c = 0; c < resiParm->nCons; c++)
        {
            CONSLISTPARMS* conslistparms = resiParm->consList[c];
            for (int s = 0; s < conslistparms->consSubListSize; s++)
            {
                CONSPARMS* consparms = conslistparms->consSubList[s];
                if (consparms->valid == 1)
                {
                    for (int i = 0; i < resiParm->nBonds; ++i)
                    {
                        BONDPARMS* bondparms = resiParm->bondList[i];

                        if (bondparms->func == 1)
                        {
                            if (consparms->atomI == bondparms->atomI && consparms->atomJ == bondparms->atomJ)
                            {
                                consparms->valid = 0;
                                break;
                            }
                            if (consparms->atomI == bondparms->atomJ && consparms->atomJ == bondparms->atomI)
                            {
                                consparms->valid = 0;
                                break;
                            }
                        }
                    }

                    for (int e = 0; e < resiParm->nExclude; e++)
                    {
                        EXCLUDEPARMS* excludeparms = resiParm->exclusionList[e];
                        if (excludeparms->valid == 1)
                        {
                            if (consparms->atomI == excludeparms->atomI && consparms->atomJ == excludeparms->atomJ)
                            {
                                consparms->valid = 0;
                                break;
                            }
                            if (consparms->atomI == excludeparms->atomJ && consparms->atomJ == excludeparms->atomI)
                            {
                                consparms->valid = 0;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

int genMartiniBondPair(CHARMM_PARMS* charmmParms, MMFF *mmff)
{

    validateExclusions(mmff);

    int bpairTotalSize = 0;

    //for(int r=0; r< charmmParms->resiConnSize; ++r){
    //    RESI_CONN* resiConn=charmmParms->resiConnList[r];
    //    bpairTotalSize+=resiConn->bondListSize;
    //}

    for (int r = 0; r < mmff->nResiParms; ++r)
    {
        RESIPARMS* resiParm = mmff->resiParms[r];
        for (int i = 0; i < resiParm->nBonds; ++i)
        {
            if (resiParm->bondList[i]->func == 1)
            {
                bpairTotalSize++;
            }
        }

        for (int e = 0; e < resiParm->nExclude; e++)
        {
            if (resiParm->exclusionList[e]->valid == 1)
            {
                bpairTotalSize++;
            }
        }

        for (int c = 0; c < resiParm->nCons; c++)
        {
            CONSLISTPARMS* conslistparms = resiParm->consList[c];
            for (int s = 0; s < conslistparms->consSubListSize; s++)
            {
                CONSPARMS* consparms = conslistparms->consSubList[s];
                if (consparms->valid == 1)
                {
                    bpairTotalSize++;
                }
            }
        }
    }

    BPAIR_CONN *bpairHeap = ddcMalloc(bpairTotalSize * sizeof (BPAIR_CONN));
    BPAIR_CONN **bpairPtrs = ddcMalloc(bpairTotalSize * sizeof (BPAIR_CONN*));
    for (int i = 0; i < bpairTotalSize; ++i)
    {
        bpairPtrs[i] = &bpairHeap[i];
    }

    bpairTotalSize = 0;

    for (int r = 0; r < mmff->nResiParms; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        RESIPARMS* resiParm = mmff->resiParms[r];
        resiConn->bpairList = &bpairPtrs[bpairTotalSize];
        int bpairSize = 0;

        for (int i = 0; i < resiParm->nBonds; ++i)
        {
            BONDPARMS* bondparms = resiParm->bondList[i];
            if (bondparms->func == 1)
            {
                bpairPtrs[bpairTotalSize]->atmI = bondparms->atomI;
                bpairPtrs[bpairTotalSize]->atmJ = bondparms->atomJ;
                char* atmTypeI = bondparms->atomTypeI;
                char* atmTypeJ = bondparms->atomTypeJ;
                double eps = 0.0;
                double sigma = 0.0;
                if (getMartiniLJpairParm(mmff, atmTypeI, atmTypeJ, &eps, &sigma) < 0)
                {
                    char message[80];
                    sprintf(message, "Martini BondPair error: Not able to find LJ parameter %s - %s \n", atmTypeI, atmTypeJ);
                    error_action(message, ERROR_IN("Martini BondPair", ABORT));
                }
                bpairPtrs[bpairTotalSize]->eps = eps;
                bpairPtrs[bpairTotalSize]->sigma = sigma;
                // Default no potential shift
                bpairPtrs[bpairTotalSize]->shift = 0.0;

                ++bpairSize;
                ++bpairTotalSize;
            }
        }

        for (int e = 0; e < resiParm->nExclude; e++)
        {
            EXCLUDEPARMS* excludeparms = resiParm->exclusionList[e];
            if (excludeparms->valid == 1)
            {
                bpairPtrs[bpairTotalSize]->atmI = excludeparms->atomI;
                bpairPtrs[bpairTotalSize]->atmJ = excludeparms->atomJ;
                char* atmTypeI = excludeparms->atomTypeI;
                char* atmTypeJ = excludeparms->atomTypeJ;
                double eps = 0;
                double sigma = 0;
                if (getMartiniLJpairParm(mmff, atmTypeI, atmTypeJ, &eps, &sigma) < 0)
                {
                    char message[80];
                    sprintf(message, "Martini BondPair error: Not able to find LJ parameter %s - %s \n", atmTypeI, atmTypeJ);
                    error_action(message, ERROR_IN("Martini BondPair", ABORT));
                }
                bpairPtrs[bpairTotalSize]->eps = eps;
                bpairPtrs[bpairTotalSize]->sigma = sigma;

                ++bpairSize;
                ++bpairTotalSize;
            }
        }

        for (int c = 0; c < resiParm->nCons; c++)
        {
            CONSLISTPARMS* conslistparms = resiParm->consList[c];
            for (int s = 0; s < conslistparms->consSubListSize; s++)
            {
                CONSPARMS* consparms = conslistparms->consSubList[s];
                if (consparms->valid == 1)
                {
                    bpairPtrs[bpairTotalSize]->atmI = consparms->atomI;
                    bpairPtrs[bpairTotalSize]->atmJ = consparms->atomJ;
                    char* atmTypeI = consparms->atomTypeI;
                    char* atmTypeJ = consparms->atomTypeJ;
                    double eps = 0;
                    double sigma = 0;
                    if (getMartiniLJpairParm(mmff, atmTypeI, atmTypeJ, &eps, &sigma) < 0)
                    {
                        char message[80];
                        sprintf(message, "Martini BondPair error: Not able to find LJ parameter %s - %s \n", atmTypeI, atmTypeJ);
                        error_action(message, ERROR_IN("Martini BondPair", ABORT));
                    }
                    bpairPtrs[bpairTotalSize]->eps = eps;
                    bpairPtrs[bpairTotalSize]->sigma = sigma;

                    ++bpairSize;
                    ++bpairTotalSize;
                }
            }
        }

        resiConn->bpairListSize = bpairSize;

    }

    return 0;
}

int findIndexInConsGroup(int atomID, int* atomList, int size)
{
    for (int i = 0; i < size; ++i)
    {
        if (atomID == atomList[i])
        {
            return i;
        }
    }
    char message[80];
    sprintf(message, "findIndexInConsGroup error: Not able to find %d\n", atomID);
    error_action(message, ERROR_IN("findIndexInConsGroup", ABORT));

    return -1;
}

void genConstraintWithZeroCONSPAIR(CHARMM_PARMS* charmmParms, MMFF *mmff)
{
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN *resiConn = charmmParms->resiConnList[r];
        RESIPARMS *resiParm = mmff->resiParms[r];

        int nAtoms = resiConn->atomListSize;
        int consIndex[nAtoms];
        int atomConsIndex[nAtoms];
        for (int i = 0; i < nAtoms; ++i)
        {
            consIndex[i] = -1; // -1 means atom has no constraints
            // -2 means constraints has already been counted.
            atomConsIndex[i] = -1;
        }

        int nCons = resiParm->nCons;
        //int numConsAtom=0;
        for (int c = 0; c < nCons; c++)
        {
            CONSLISTPARMS* conslistparms = resiParm->consList[c];
            //numConsAtom+=conslistparms->consSubListSize;
            for (int s = 0; s < conslistparms->consSubListSize; s++)
            {
                CONSPARMS* consparms = conslistparms->consSubList[s];
                consIndex[consparms->atomI] = c;
                consIndex[consparms->atomJ] = c;
                atomConsIndex[consparms->atomI] = c;
                atomConsIndex[consparms->atomJ] = c;
            }
        }

        int consGroup[nCons];
        int numAtomInGroup[nCons];
        for (int i = 0; i < nCons; ++i)
        {
            consGroup[i] = -1; // -1 the group has not appeared in the consIndex yet
            numAtomInGroup[i] = 0; // For count the number of atoms in the constraint group
        }

        for (int i = 0; i < nAtoms; ++i)
        {
            if (consIndex[i] >= 0)
            {
                numAtomInGroup[consIndex[i]]++;
            }
        }

        //int totalConsAtom=0;
        //for(int i=0; i<nCons; ++i){
        //    totalConsAtom+=numAtomInGroup[i];
        //}

        for (int i = 0; i < nAtoms; ++i)
        {
            if (consIndex[i] >= 0)
            {
                if (consGroup[consIndex[i]] > 0)
                {
                    consIndex[i] = -2; // constraints has already been counted label -2 for skipping it
                }
                else
                {
                    consGroup[consIndex[i]] = 1;
                }
            }
        }

        int nConnCons = 0;
        //int testConsGroup=0;
        for (int i = 0; i < nAtoms; ++i)
        {
            if (consIndex[i] >= -1)
            {
                nConnCons++;
            }
            //if(consIndex[i]>=0){
            //    testConsGroup++;
            //}             
        }

        CONSTRAINT* connConsStack = ddcMalloc(nConnCons * sizeof (CONSTRAINT));
        resiConn->consList = ddcMalloc(nConnCons * sizeof (CONSTRAINT*));
        resiConn->consListSize = nConnCons;
        for (int i = 0; i < nConnCons; ++i)
        {
            resiConn->consList[i] = &connConsStack[i];
        }

        nConnCons = 0;
        for (int i = 0; i < nAtoms; ++i)
        {
            int consIndexCur = consIndex[i];
            if (consIndexCur >= 0)
            {
                CONSTRAINT* cons = resiConn->consList[nConnCons];
                TATOM_PARMS* atom = resiConn->atomList[i];
                cons->atomID = atom->atmID;

                CONSLISTPARMS* conslistparms = resiParm->consList[consIndexCur];
                cons->numPair = conslistparms->consSubListSize;
                cons->numAtom = numAtomInGroup[consIndex[i]];
                CONS_PAIR* consPairStack = ddcMalloc(cons->numPair * sizeof (CONS_PAIR));
                cons->conspairList = ddcMalloc(cons->numPair * sizeof (CONS_PAIR*));
                cons->atomIDList = ddcMalloc(cons->numAtom * sizeof (int));
                int atomIDcount = 0;
                for (int a = 0; a < nAtoms; ++a)
                {
                    if (atomConsIndex[a] == consIndexCur)
                    {
                        assert(atomIDcount < cons->numAtom); //Make sure wont go over the boundary
                        cons->atomIDList[atomIDcount] = a;
                        ++atomIDcount;
                    }
                }
                for (int c = 0; c < cons->numPair; ++c)
                {
                    cons->conspairList[c] = &consPairStack[c];
                    int atomI = conslistparms->consSubList[c]->atomI;
                    cons->conspairList[c]->atomI = atomI;
                    cons->conspairList[c]->atomIindex = findIndexInConsGroup(atomI, cons->atomIDList, cons->numAtom);
                    int atomJ = conslistparms->consSubList[c]->atomJ;
                    cons->conspairList[c]->atomJ = atomJ;
                    cons->conspairList[c]->atomJindex = findIndexInConsGroup(atomJ, cons->atomIDList, cons->numAtom);
                    cons->conspairList[c]->distance = conslistparms->consSubList[c]->r0;
                }

                nConnCons++;
            }
            else if (consIndexCur == -1)
            {
                CONSTRAINT* cons = resiConn->consList[nConnCons];
                TATOM_PARMS* atom = resiConn->atomList[i];
                cons->atomID = atom->atmID;
                cons->numPair = 0;
                cons->numAtom = 0;
                cons->conspairList = NULL;
                nConnCons++;
            }
        }

    }
}

void genConstraint(CHARMM_PARMS* charmmParms, MMFF *mmff)
{
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN *resiConn = charmmParms->resiConnList[r];
        RESIPARMS *resiParm = mmff->resiParms[r];

        int nAtoms = resiConn->atomListSize;
        int consIndex[nAtoms];
        int atomConsIndex[nAtoms];
        for (int i = 0; i < nAtoms; ++i)
        {
            consIndex[i] = -1; // -1 means atom has no constraints
            // -2 means constraints has already been counted.
            atomConsIndex[i] = -1;
        }

        int nCons = resiParm->nCons;
        //int numConsAtom=0;
        for (int c = 0; c < nCons; c++)
        {
            CONSLISTPARMS* conslistparms = resiParm->consList[c];
            //numConsAtom+=conslistparms->consSubListSize;
            for (int s = 0; s < conslistparms->consSubListSize; s++)
            {
                CONSPARMS* consparms = conslistparms->consSubList[s];
                consIndex[consparms->atomI] = c;
                consIndex[consparms->atomJ] = c;
                atomConsIndex[consparms->atomI] = c;
                atomConsIndex[consparms->atomJ] = c;
            }
        }

        int consGroup[nCons];
        int numAtomInGroup[nCons];
        for (int i = 0; i < nCons; ++i)
        {
            consGroup[i] = -1; // -1 the group has not appeared in the consIndex yet
            numAtomInGroup[i] = 0; // For count the number of atoms in the constraint group
        }

        for (int i = 0; i < nAtoms; ++i)
        {
            if (consIndex[i] >= 0)
            {
                numAtomInGroup[consIndex[i]]++;
            }
        }

        //int totalConsAtom=0;
        //for(int i=0; i<nCons; ++i){
        //    totalConsAtom+=numAtomInGroup[i];
        //}

        for (int i = 0; i < nAtoms; ++i)
        {
            if (consIndex[i] >= 0)
            {
                if (consGroup[consIndex[i]] > 0)
                {
                    consIndex[i] = -2; // constraints has already been counted label -2 for skipping it
                }
                else
                {
                    consGroup[consIndex[i]] = 1;
                }
            }
        }

        CONSTRAINT* connConsStack = ddcMalloc(nCons * sizeof (CONSTRAINT));
        resiConn->consList = ddcMalloc(nCons * sizeof (CONSTRAINT*));
        resiConn->consListSize = nCons;
        for (int i = 0; i < nCons; ++i)
        {
            resiConn->consList[i] = &connConsStack[i];
        }

        int nConnCons = 0;
        for (int i = 0; i < nAtoms; ++i)
        {
            int consIndexCur = consIndex[i];
            if (consIndexCur >= 0)
            {
                CONSTRAINT* cons = resiConn->consList[nConnCons];
                TATOM_PARMS* atom = resiConn->atomList[i];
                cons->atomID = atom->atmID;

                CONSLISTPARMS* conslistparms = resiParm->consList[consIndexCur];
                cons->numPair = conslistparms->consSubListSize;
                cons->numAtom = numAtomInGroup[consIndex[i]];
                CONS_PAIR* consPairStack = ddcMalloc(cons->numPair * sizeof (CONS_PAIR));
                cons->conspairList = ddcMalloc(cons->numPair * sizeof (CONS_PAIR*));
                cons->atomIDList = ddcMalloc(cons->numAtom * sizeof (int));
                int atomIDcount = 0;
                for (int a = 0; a < nAtoms; ++a)
                {
                    if (atomConsIndex[a] == consIndexCur)
                    {
                        assert(atomIDcount < cons->numAtom); //Make sure wont go over the boundary
                        cons->atomIDList[atomIDcount] = a;
                        ++atomIDcount;
                    }
                }
                for (int c = 0; c < cons->numPair; ++c)
                {
                    cons->conspairList[c] = &consPairStack[c];
                    int atomI = conslistparms->consSubList[c]->atomI;
                    cons->conspairList[c]->atomI = atomI;
                    cons->conspairList[c]->atomIindex = findIndexInConsGroup(atomI, cons->atomIDList, cons->numAtom);
                    int atomJ = conslistparms->consSubList[c]->atomJ;
                    cons->conspairList[c]->atomJ = atomJ;
                    cons->conspairList[c]->atomJindex = findIndexInConsGroup(atomJ, cons->atomIDList, cons->numAtom);
                    cons->conspairList[c]->distance = conslistparms->consSubList[c]->r0;
                }

                nConnCons++;
            }
        }

    }
}

int genMartiniConn(CHARMM_PARMS *charmmParms, MMFF *mmff)
{
    charmmParms->resiConnSize = mmff->nResiParms;

    // MASSPARMS
    charmmParms->charmmTop = ddcMalloc(sizeof (CHARMM_TOP));
    MASS_PARMS* massStack = ddcMalloc(mmff->nAtomType * sizeof (MASS_PARMS));
    charmmParms->charmmTop->massParms = ddcMalloc(mmff->nAtomType * sizeof (MASS_PARMS*));
    charmmParms->charmmTop->massParmSize = mmff->nAtomType;
    for (int i = 0; i < mmff->nAtomType; ++i)
    {
        charmmParms->charmmTop->massParms[i] = &massStack[i];
        strcpy(charmmParms->charmmTop->massParms[i]->atmType, mmff->atomTypes[i]->atomType);
        charmmParms->charmmTop->massParms[i]->atmTypeID = mmff->atomTypes[i]->atomTypeID;
        charmmParms->charmmTop->massParms[i]->mass = mmff->atomTypes[i]->mass;
    }

    RESI_CONN *resConnStack = ddcMalloc(charmmParms->resiConnSize * sizeof (RESI_CONN));
    charmmParms->resiConnList = ddcMalloc(charmmParms->resiConnSize * sizeof (RESI_CONN*));

    for (int i = 0; i < charmmParms->resiConnSize; ++i)
    {
        charmmParms->resiConnList[i] = &resConnStack[i];
        RESI_CONN *resiConn = charmmParms->resiConnList[i];
        RESIPARMS *resiParm = mmff->resiParms[i];
        resiConn->resID = resiParm->resID;
        resiConn->resType = resiParm->resType;
        resiConn->charge = resiParm->charge;
        resiConn->centerAtom = resiParm->centerAtom;

        strcpy(resiConn->resName, resiParm->resName);
        resiConn->nTer = 0;
        resiConn->cTer = 0;

        // Group & Atom
        resiConn->groupListSize = resiParm->nGroups;
        TGROUP_PARMS* groupStack = ddcMalloc(resiConn->groupListSize * sizeof (TGROUP_PARMS));
        resiConn->groupList = ddcMalloc(resiConn->groupListSize * sizeof (TGROUP_PARMS*));
        int totResAtom = 0;
        for (int j = 0; j < resiConn->groupListSize; j++)
        {
            totResAtom = totResAtom + resiParm->groupList[j]->nAtoms;
        }
        // Number of atom should be always large than 0;
        assert(totResAtom > 0);
        resiConn->atomListSize = totResAtom;
        TATOM_PARMS* atomStack = ddcMalloc(totResAtom * sizeof (TATOM_PARMS));
        resiConn->atomList = ddcMalloc(totResAtom * sizeof (TATOM_PARMS*));

        //for (int j = 0; j < totResAtom; j++) {
        //    resiConn->atomList[j] = &atomStack[j];
        //}

        int index = 0;
        for (int j = 0; j < resiConn->groupListSize; j++)
        {
            resiConn->groupList[j] = &groupStack[j];
            resiConn->groupList[j]->grpID = resiParm->groupList[j]->groupID;
            resiConn->groupList[j]->grpAtomSize = resiParm->groupList[j]->nAtoms;
            resiConn->groupList[j]->grpAtoms = ddcMalloc(resiParm->groupList[j]->nAtoms * sizeof (TATOM_PARMS*));

            for (int k = 0; k < resiParm->groupList[j]->nAtoms; k++)
            {
                resiConn->atomList[index] = &atomStack[index];
                resiConn->atomList[index]->atmID = resiParm->groupList[j]->atomList[k]->atomID;
                strcpy(resiConn->atomList[index]->atmName, resiParm->groupList[j]->atomList[k]->atomName);
                resiConn->atomList[index]->charge = resiParm->groupList[j]->atomList[k]->charge;
                resiConn->atomList[index]->atmTypePtr = charmmParms->charmmTop->massParms[resiParm->groupList[j]->atomList[k]->atomTypeID];

                resiConn->groupList[j]->grpAtoms[k] = resiConn->atomList[index];

                index++;
            }
        }

        // BOND
        resiConn->bondListSize = resiParm->nBonds;
        if (resiConn->bondListSize > 0)
        {
            BOND_CONN* bondConnStack = ddcMalloc(resiConn->bondListSize * sizeof (BOND_CONN));
            resiConn->bondList = ddcMalloc(resiConn->bondListSize * sizeof (BOND_CONN*));
            BOND_PARMS* bondParmStack = ddcMalloc(resiConn->bondListSize * sizeof (BOND_PARMS));
            for (int j = 0; j < resiConn->bondListSize; j++)
            {
                resiConn->bondList[j] = &bondConnStack[j];
                resiConn->bondList[j]->atmI = resiParm->bondList[j]->atomI;
                resiConn->bondList[j]->atmJ = resiParm->bondList[j]->atomJ;
                resiConn->bondList[j]->bondPtr = &bondParmStack[j];
                resiConn->bondList[j]->bondPtr->kb = resiParm->bondList[j]->kb;
                resiConn->bondList[j]->bondPtr->b0 = resiParm->bondList[j]->b0;
                strcpy(resiConn->bondList[j]->bondPtr->atmTypeI, resiParm->bondList[j]->atomTypeI);
                strcpy(resiConn->bondList[j]->bondPtr->atmTypeJ, resiParm->bondList[j]->atomTypeJ);
            }
        }

        // ANGLE, Cosine-based ANGLE, restrained bending potential (resanlge)
        int angleListSize = 0;
        int cosangleListSize = 0;
        int rebangleListSize = 0;
        for (int j = 0; j < resiParm->nAngles; j++)
        {
            if (resiParm->angleList[j]->func == 1)
            {
                angleListSize++;
            }
            if (resiParm->angleList[j]->func == 2)
            {
                cosangleListSize++;
            }
            if (resiParm->angleList[j]->func == 10)
            {
                rebangleListSize++;
            }            
        }
        resiConn->angleListSize = angleListSize;
        if (resiConn->angleListSize > 0)
        {
            ANGLE_CONN* angleConnStack = ddcMalloc(angleListSize * sizeof (ANGLE_CONN));
            resiConn->angleList = ddcMalloc(angleListSize * sizeof (ANGLE_CONN*));
            ANGLE_PARMS* angleParmStack = ddcMalloc(angleListSize * sizeof (ANGLE_PARMS));
            int ind = 0;
            for (int j = 0; j < resiParm->nAngles; j++)
            {
                if (resiParm->angleList[j]->func == 1)
                {
                    resiConn->angleList[ind] = &angleConnStack[ind];
                    resiConn->angleList[ind]->atmI = resiParm->angleList[j]->atomI;
                    resiConn->angleList[ind]->atmJ = resiParm->angleList[j]->atomJ;
                    resiConn->angleList[ind]->atmK = resiParm->angleList[j]->atomK;
                    resiConn->angleList[ind]->anglePtr = &angleParmStack[ind];
                    resiConn->angleList[ind]->anglePtr->ktheta = resiParm->angleList[j]->ktheta;
                    resiConn->angleList[ind]->anglePtr->theta0 = resiParm->angleList[j]->theta0;
                    resiConn->angleList[ind]->anglePtr->kub = 0;
                    resiConn->angleList[ind]->anglePtr->s0 = 0;
                    ind++;
                }
            }
        }

        resiConn->cosangleListSize = cosangleListSize;
        if (resiConn->cosangleListSize > 0)
        {
            ANGLE_CONN* cosangleConnStack = ddcMalloc(cosangleListSize * sizeof (ANGLE_CONN));
            resiConn->cosangleList = ddcMalloc(cosangleListSize * sizeof (ANGLE_CONN*));
            ANGLE_PARMS* cosangleParmStack = ddcMalloc(cosangleListSize * sizeof (ANGLE_PARMS));
            int ind = 0;
            for (int j = 0; j < resiParm->nAngles; j++)
            {
                if (resiParm->angleList[j]->func == 2)
                {
                    resiConn->cosangleList[ind] = &cosangleConnStack[ind];
                    resiConn->cosangleList[ind]->atmI = resiParm->angleList[j]->atomI;
                    resiConn->cosangleList[ind]->atmJ = resiParm->angleList[j]->atomJ;
                    resiConn->cosangleList[ind]->atmK = resiParm->angleList[j]->atomK;
                    resiConn->cosangleList[ind]->anglePtr = &cosangleParmStack[ind];
                    resiConn->cosangleList[ind]->anglePtr->ktheta = resiParm->angleList[j]->ktheta;
                    resiConn->cosangleList[ind]->anglePtr->theta0 = resiParm->angleList[j]->theta0;
                    resiConn->cosangleList[ind]->anglePtr->kub = 0;
                    resiConn->cosangleList[ind]->anglePtr->s0 = 0;
                    ind++;
                }
            }
        }

        resiConn->rebangleListSize = rebangleListSize;
        if (resiConn->rebangleListSize > 0)
        {
            ANGLE_CONN* rebangleConnStack = ddcMalloc(rebangleListSize * sizeof (ANGLE_CONN));
            resiConn->rebangleList = ddcMalloc(rebangleListSize * sizeof (ANGLE_CONN*));
            ANGLE_PARMS* rebangleParmStack = ddcMalloc(rebangleListSize * sizeof (ANGLE_PARMS));
            int ind = 0;
            for (int j = 0; j < resiParm->nAngles; j++)
            {
                if (resiParm->angleList[j]->func == 10)
                {
                    resiConn->rebangleList[ind] = &rebangleConnStack[ind];
                    resiConn->rebangleList[ind]->atmI = resiParm->angleList[j]->atomI;
                    resiConn->rebangleList[ind]->atmJ = resiParm->angleList[j]->atomJ;
                    resiConn->rebangleList[ind]->atmK = resiParm->angleList[j]->atomK;
                    resiConn->rebangleList[ind]->anglePtr = &rebangleParmStack[ind];
                    resiConn->rebangleList[ind]->anglePtr->ktheta = resiParm->angleList[j]->ktheta;
                    resiConn->rebangleList[ind]->anglePtr->theta0 = resiParm->angleList[j]->theta0;
                    resiConn->rebangleList[ind]->anglePtr->kub = 0;
                    resiConn->rebangleList[ind]->anglePtr->s0 = 0;
                    ind++;
                }
            }
        }
        
        // DIHEDRAL TORSION and IMPROPER
        int torsListSize = 0;
        int imprListSize = 0;
        for (int j = 0; j < resiParm->nTors; j++)
        {
            if (resiParm->torsList[j]->func == 1)
            {
                torsListSize++;
            }
            if (resiParm->torsList[j]->func == 2)
            {
                imprListSize++;
            }
        }
        resiConn->torsListSize = torsListSize;
        if (resiConn->torsListSize > 0)
        {
            TORS_CONN* torsConnStack = ddcMalloc(torsListSize * sizeof (TORS_CONN));
            resiConn->torsList = ddcMalloc(torsListSize * sizeof (TORS_CONN*));
            TORSION_PARMS* torsParmStack = ddcMalloc(torsListSize * sizeof (TORSION_PARMS));
            int ind = 0;
            for (int j = 0; j < resiParm->nTors; j++)
            {
                if (resiParm->torsList[j]->func == 1)
                {
                    resiConn->torsList[ind] = &torsConnStack[ind];
                    resiConn->torsList[ind]->atmI = resiParm->torsList[j]->atomI;
                    resiConn->torsList[ind]->atmJ = resiParm->torsList[j]->atomJ;
                    resiConn->torsList[ind]->atmK = resiParm->torsList[j]->atomK;
                    resiConn->torsList[ind]->atmL = resiParm->torsList[j]->atomL;
                    resiConn->torsList[ind]->torsPtr = &torsParmStack[ind];
                    resiConn->torsList[ind]->torsPtr->kchi = resiParm->torsList[j]->kchi;
                    resiConn->torsList[ind]->torsPtr->delta = resiParm->torsList[j]->delta;
                    resiConn->torsList[ind]->torsPtr->n = resiParm->torsList[j]->n;
                    ind++;
                }
            }
        }

        resiConn->imprListSize = imprListSize;
        if (resiConn->imprListSize > 0)
        {
            IMPR_CONN* imprConnStack = ddcMalloc(imprListSize * sizeof (IMPR_CONN));
            resiConn->imprList = ddcMalloc(imprListSize * sizeof (IMPR_CONN*));
            IMPROPER_PARMS* imprParmStack = ddcMalloc(imprListSize * sizeof (IMPROPER_PARMS));
            int ind = 0;
            for (int j = 0; j < resiParm->nTors; j++)
            {
                if (resiParm->torsList[j]->func == 2)
                {
                    resiConn->imprList[ind] = &imprConnStack[ind];
                    resiConn->imprList[ind]->atmI = resiParm->torsList[j]->atomI;
                    resiConn->imprList[ind]->atmJ = resiParm->torsList[j]->atomJ;
                    resiConn->imprList[ind]->atmK = resiParm->torsList[j]->atomK;
                    resiConn->imprList[ind]->atmL = resiParm->torsList[j]->atomL;
                    resiConn->imprList[ind]->imprPtr = &imprParmStack[ind];
                    resiConn->imprList[ind]->imprPtr->kpsi = resiParm->torsList[j]->kchi;
                    resiConn->imprList[ind]->imprPtr->psi0 = resiParm->torsList[j]->delta;
                    ind++;
                }
            }
        }

        // The rest lists are set to 0
        resiConn->cmapListSize = 0;
        resiConn->cmapList = NULL;
        resiConn->weightListSize = 0;
        resiConn->weightList = NULL;

    }

    sortBondList(charmmParms);

    // BondPair
    genMartiniBondPair(charmmParms, mmff);

    genAtmRange(charmmParms);

    // Generate constraint
    genConstraint(charmmParms, mmff);

    return 0;
}

double CGLennardJones_setShift(double sigma, double eps, double rcut)
{
    double sigma_r = sigma / rcut;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;
    return (-4.0 * eps * (s12 - s6));
}

void BpairLennardJones_setShift(CHARMMPOT_PARMS *parms)
{
    POTENTIAL *potential = parms->potential;
    double cutoff;
    object_get((OBJECT*) potential, "cutoff", &cutoff, WITH_UNITS, 1, "11.0", "Angstrom", NULL);
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        for (int b = 0; b < resiConn->bpairListSize; ++b)
        {
            BPAIR_CONN* bpairConn = resiConn->bpairList[b];
            bpairConn->shift = CGLennardJones_setShift(bpairConn->sigma, bpairConn->eps, cutoff);
        }
    }

}

CharmmLJ_PARMS **martiniLJ_parms(CHARMMPOT_PARMS *parms, MMFF* mmff)
{
    POTENTIAL *potential = parms->potential;
    int nspecies = parms->nspecies;
    //parms->charmmParms->charmmPar->ljParmSize = nspecies;
    double cutoff;
    int potentialShift;
    CharmmLJ_PARMS **ljParms = (CharmmLJ_PARMS **) ddcMalloc(sizeof (CharmmLJ_PARMS*) * nspecies * nspecies);
    object_get((OBJECT*) potential, "cutoff", &cutoff, WITH_UNITS, 1, "11.0", "Angstrom", NULL);
    object_get((OBJECT*) potential, "potential-shift", &potentialShift, INT, 1, "1");
    //parms->nonBond_fcn = (double(*)(void *, double, double *))CHLennardJones;
    parms->rmax = cutoff;

    if (potentialShift)
    {
        BpairLennardJones_setShift(parms);
    }



    int nlist = nspecies + nspecies * (nspecies - 1) / 2;

    CharmmLJ_PARMS *buffer = ddcMalloc(nlist * sizeof (CharmmLJ_PARMS));

    // buffer matrix looks like:
    //      0   1   2   3   4
    //  0   0   5   9   12  14   
    //  1       1   6   10  13
    //  2           2   7   11
    //  3               3   8
    //  4                   4
    for (int i = 0; i < mmff->nLJParms; i++)
    {
        LJPARMS* mmff_lj = mmff->ljParms[i];
        CharmmLJ_PARMS *ljParms_ab;
        int i;
        int j;
        if (mmff_lj->indexI > mmff_lj->indexJ)
        {
            i = mmff_lj->indexI;
            j = mmff_lj->indexJ;
        }
        else
        {
            i = mmff_lj->indexJ;
            j = mmff_lj->indexI;
        }
        int d = i - j;
        int offset = 0;
        for (int i = 0; i < d; i++)
        {
            offset += nspecies - i;
        }

        //if (mmff_lj->indexI == mmff_lj->indexJ) {
        //    ljParms_ab = buffer + mmff_lj->indexI;
        //} else if(mmff_lj->indexI > mmff_lj->indexJ) {
        //    ljParms_ab = buffer + mmff_lj->indexI + mmff_lj->indexJ*nspecies;
        //} else{
        //    ljParms_ab = buffer + mmff_lj->indexJ + mmff_lj->indexI*nspecies;
        //}
        ljParms_ab = buffer + offset + j;

        ljParms_ab->eps = mmff_lj->eps;
        ljParms_ab->sigma = mmff_lj->sigma;
        ljParms_ab->rcut = cutoff;
        if (potentialShift)
        {
            ljParms_ab->shift = CGLennardJones_setShift(ljParms_ab->sigma, ljParms_ab->eps, ljParms_ab->rcut);
        }
        else
        {
            ljParms_ab->shift = 0.0;
        }

        int ab = mmff_lj->indexI + mmff_lj->indexJ*nspecies;
        int ba = mmff_lj->indexJ + mmff_lj->indexI*nspecies;

        ljParms[ab] = ljParms[ba] = ljParms_ab;
    }

    return ljParms;
}

int getCGLJindexbySpecie(SPECIES* specie, CHARMM_PARMS *charmmParms)
{


    //int index=0;
    char resName[5];

    char* specieStr = specie->name;

    char* delimPtr = strpbrk(specieStr, "cnx");
    //char terminus = *delimPtr;
    int resNSize = delimPtr - specieStr;
    strncpy(resName, specieStr, resNSize);
    resName[resNSize] = '\0';
    char *atmName = delimPtr + 1;

    RESI_CONN** resConnList = charmmParms->resiConnList;
    int resiConnSize = charmmParms->resiConnSize;
    for (int i = 0; i < resiConnSize; ++i)
    {
        if (resNameDiffer(resConnList[i]->resName, resName) == 0)
        {
            TATOM_PARMS** atmList = resConnList[i]->atomList;
            int atmListSize = resConnList[i]->atomListSize;
            for (int j = 0; j < atmListSize; ++j)
            {
                if (strcmp(atmList[j]->atmName, atmName) == 0)
                {
                    return atmList[j]->atmTypePtr->atmTypeID;
                }
            }
        }
    }
    error_action("getCGLJindexbySpecie error: Not able to find index by specie", ERROR_IN("getCGLJindexbySpecie", ABORT));
    return -1;
}

void martiniNonBond(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    profile(CHARMM_NONBOND, START);
    //int IndexRmax = -1;
    //for (int i = 0; i < sys->neighbor->nc; i++) {
    //   if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
    //}
    //if (IndexRmax == -1) 
    //  error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));

    CHARMM_PARMS *charmmParms = parms->charmmParms;
    int nspecies = parms->nspecies; // has been set to mmff->nAtomType
    NPARTICLE *particles = sys->neighbor->particles;
    double r2cut = SQ(parms->rmax);
    //double rcoulomb=parms->rcoulomb;
    double krf = parms->krf;
    double crf = parms->crf;

    STATE *state = sys->collection->state;
    double *rx = state->rx;
    double *ry = state->ry;
    double *rz = state->rz;
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;
    double *q = state->q;
    double minspan;
    box_get(NULL, MINSPAN, &minspan);
    double R2cut = 0.25 * minspan*minspan;

    int speciesIndex[sys->nspecies];
    for (int j = 0; j < sys->nspecies; j++)
    {
        speciesIndex[sys->species[j]->index] = getCGLJindexbySpecie(sys->species[j], charmmParms);
    }
    int sIndex[sys->nion];
    THREE_SMATRIX virial = szero;
    for (unsigned i = 0; i < sys->nion; i++)
    {
        sIndex[i] = speciesIndex[state->species[i]->index];
    }
    double q2 = 0.0;
    for (unsigned i = 0; i < sys->nlocal; i++) q2 += SQ(q[i]); // Should be in the sys->local or sys->nion?

    double keR = ke / parms->epsilon_r;
    double vLJ = 0.0;
    double vEle = -0.5 * q2 * keR*crf; //Self electronic energy
    for (unsigned i = 0; i < sys->nlocal; i++)
    {
        int si = sIndex[i];
        double kqi = keR * q[i];
        double xi = rx[i];
        double yi = ry[i];
        double zi = rz[i];
        double fxi = 0.0;
        double fyi = 0.0;
        double fzi = 0.0;

        PAIRS *pij = particles[i].ifirst[0];
        while (pij != NULL)
        {
            int j = pij->j;
            int sj = sIndex[j];
            int sij = sj + nspecies*si;
            //double rcut = parms->rcut[sij].value;
            double x = xi - rx[j];
            double y = yi - ry[j];
            double z = zi - rz[j];
            double r2 = x * x + y * y + z*z;
            if (r2 > R2cut)
            {
                nearestImage_fast(&x, &y, &z);
                r2 = x * x + y * y + z*z;
            }
            if (r2 < r2cut)
            {
                CharmmLJ_PARMS *ljparms = parms->parms[sij];
                double sigma = ljparms->sigma;
                double eps = ljparms->eps;

                double ir = sqrt(1.0 / r2);
                double ir2 = ir*ir;
                //double r = r2*ir; 

                double sigma_r = sigma * ir;
                double s2 = sigma_r*sigma_r;
                double s4 = s2*s2;
                double s6 = s4*s2;
                double s12 = s6*s6;
                vLJ += 4.0 * eps * (s12 - s6) + ljparms->shift;
                //double vLJij =   4.0 * eps * (s12 - s6) + ljparms->shift;
                double dvdr = 24.0 * eps * (s6 - 2.0 * s12) * ir2;

                double kqij = kqi * q[j];
                vEle += kqij * (ir + krf * r2 - crf);
                //double vEleij= kqij * (ir + krf*r2 - crf);
                dvdr += kqij * (2 * krf - ir2 * ir); // F=-ke*qi*qj/er*[1/r^2 - 2krf*r]r^/r

                double fxij = -dvdr * x;
                double fyij = -dvdr * y;
                double fzij = -dvdr * z;

                fxi += fxij;
                fyi += fyij;
                fzi += fzij;
                fx[j] -= fxij;
                fy[j] -= fyij;
                fz[j] -= fzij;

                virial.xx += fxij*x;
                virial.yy += fyij*y;
                virial.zz += fzij*z;
                virial.xy += fxij*y;
                virial.xz += fxij*z;
                virial.yz += fyij*z;
            }
            pij = pij->ilink;
        }
        fx[i] += fxi;
        fy[i] += fyi;
        fz[i] += fzi;
    }
    parms->bioEnergies.nonBond = (vLJ + vEle);
    parms->bioEnergies.lj = vLJ;
    parms->bioEnergies.ele = vEle;
    e->eion += (vLJ + vEle);
    e->virial.xx += virial.xx;
    e->virial.yy += virial.yy;
    e->virial.zz += virial.zz;
    e->virial.xy += virial.xy;
    e->virial.xz += virial.xz;
    e->virial.yz += virial.yz;
    profile(CHARMM_NONBOND, END);
}

void martiniIntraMoleReaction(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    //CHARMM_PARMS *charmmParms=parms->charmmParms;
    NPARTICLE *particles = sys->neighbor->particles;
    double r2cut = SQ(parms->rmax);
    double krf = parms->krf;
    double crf = parms->crf;

    STATE *state = sys->collection->state;
    double *rx = state->rx;
    double *ry = state->ry;
    double *rz = state->rz;
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;
    double *q = state->q;
    double minspan;
    box_get(NULL, MINSPAN, &minspan);
    double R2cut = 0.25 * minspan*minspan;

    THREE_SMATRIX virial = szero;
    double keR = ke / parms->epsilon_r;
    double vEle = 0.0;
    for (unsigned i = 0; i < sys->nlocal; i++)
    {
        double kqi = keR * q[i];
        double xi = rx[i];
        double yi = ry[i];
        double zi = rz[i];
        double fxi = 0.0;
        double fyi = 0.0;
        double fzi = 0.0;

        PAIRS *pij = particles[i].ifirst[1];
        while (pij != NULL)
        {
            int j = pij->j;
            double x = xi - rx[j];
            double y = yi - ry[j];
            double z = zi - rz[j];
            double r2 = x * x + y * y + z*z;
            if (r2 > R2cut)
            {
                nearestImage_fast(&x, &y, &z);
                r2 = x * x + y * y + z*z;
            }
            if (r2 < r2cut)
            {
                double kqij = kqi * q[j];
                vEle += kqij * (krf * r2 - crf);
                double dvdr = kqij * (2 * krf);
                double fxij = -dvdr * x;
                double fyij = -dvdr * y;
                double fzij = -dvdr * z;

                fxi += fxij;
                fyi += fyij;
                fzi += fzij;
                fx[j] -= fxij;
                fy[j] -= fyij;
                fz[j] -= fzij;

                virial.xx += fxij*x;
                virial.yy += fyij*y;
                virial.zz += fzij*z;
                virial.xy += fxij*y;
                virial.xz += fxij*z;
                virial.yz += fyij*z;
            }
            pij = pij->ilink;
        }
        fx[i] += fxi;
        fy[i] += fyi;
        fz[i] += fzi;
    }
    parms->bioEnergies.nonBond += vEle;
    parms->bioEnergies.ele += vEle;
    e->eion += vEle;
    e->virial.xx += virial.xx;
    e->virial.yy += virial.yy;
    e->virial.zz += virial.zz;
    e->virial.xy += virial.xy;
    e->virial.xz += virial.xz;
    e->virial.yz += virial.yz;
}

CHARMMPOT_PARMS *martini_parms(POTENTIAL *potential)
{
    timestamp("Start Bio Martini Setup");
    CHARMMPOT_PARMS* parms = ddcMalloc(sizeof (CHARMMPOT_PARMS));
    parms->charmmParms = ddcMalloc(sizeof (CHARMM_PARMS));
    char *parmfile;
    object_get((OBJECT *) potential, "parmfile", &parmfile, STRING, 1, "martini.data");
    object_get((OBJECT *) potential, "excludePotentialTerm", &parms->excludePotentialTerm, INT, 1, "0");

    object_compilefile(parmfile);

    MMFF *mmff = mmff_init(potential, "martini");

    genMartiniConn(parms->charmmParms, mmff);

    parms->nspecies = mmff->nAtomType;
    parms->potential = potential;

    parms->parms = (void *) martiniLJ_parms(parms, mmff);

    object_get((OBJECT*) potential, "rmax4all", &(parms->rmax4all), WITH_UNITS, 1, "11.0", "Angstrom", NULL); // Martini use 11 angstrom cutoff by default
    object_get((OBJECT *) potential, "rcoulomb", &(parms->rcoulomb), WITH_UNITS, 1, "11.0", "Angstrom", NULL);
    object_get((OBJECT *) potential, "epsilon_r", &(parms->epsilon_r), DOUBLE, 1, "15.0");
    object_get((OBJECT *) potential, "epsilon_rf", &(parms->epsilon_rf), DOUBLE, 1, "-1.0");
    double irc = 1.0 / parms->rcoulomb;
    double irc3 = irc * irc*irc;
    if (parms->epsilon_rf != -1.0)
    {
        parms->krf = (parms->epsilon_rf - parms->epsilon_r) / (2 * parms->epsilon_rf + parms->epsilon_r) * irc3;
        parms->crf = 3 * (parms->epsilon_rf) / (2 * parms->epsilon_rf + parms->epsilon_r) * irc;
    }
    else
    {
        parms->krf = 0.5 * irc3;
        parms->crf = 1.5 * irc;
    }

    parms->gidOrder = NULL;
    parms->residueSet.list = NULL;
    parms->residueSet.molSize = 0;
    parms->residueSet.listSize = 0;
    parms->statechpad.rx = NULL;
    parms->statechpad.ry = NULL;
    parms->statechpad.rz = NULL;
    parms->statechpad.vx = NULL;
    parms->statechpad.vy = NULL;
    parms->statechpad.vz = NULL;
    parms->statechpad.fx = NULL;
    parms->statechpad.fy = NULL;
    parms->statechpad.fz = NULL;
    parms->statechpad.label = NULL;
    parms->statechpad.species = NULL;

    if (parms->rmax > parms->rmax4all) parms->rmax4all = parms->rmax;
    if (parms->rcoulomb > parms->rmax4all) parms->rmax4all = parms->rcoulomb;

    SYSTEM *sys = NULL;
    sys = system_getSystem(sys);
    //look for species that determine task->Residue ownership
    int isResidueOwner = 0;
    int *isResidueOwnerArr = (int*) ddcMalloc(sizeof (int)*sys->nspecies);
    object_get((OBJECT *) potential, "useAutoResidueOwnership", &(parms->useAutoResidueOwnership), INT, 1, "1");
    char speciesName[64];
    parms->speciesInfo = (SPECIESINFO *) ddcMalloc(sys->nspecies * sizeof (SPECIESINFO));
    for (int i = 0; i < sys->nspecies; i++)
    {
        char * spname = sys->species[i]->name;
        OBJECT * spTemp = object_find(spname, "SPECIES");
        object_get((OBJECT *) spTemp, "isResidueOwner", &isResidueOwner, INT, 1, "0");
        isResidueOwnerArr[i] = isResidueOwner;
        strcpy(speciesName, spname);
        char *delim = index(speciesName, 'x');
        *delim = (char) 0;
        char *residueName = speciesName;
        char *atomName = delim + 1;
        int nResi = mmff->nResiParms;
        int rIndex = 0;
        for (; rIndex < nResi; rIndex++) if (strcmp(residueName, mmff->resiParms[rIndex]->name) == 0) break;
        assert(rIndex != nResi);
        RESIPARMS *resiParms = mmff->resiParms[rIndex];
        int offset = 0;
        for (int gIndex = 0; gIndex < resiParms->nGroups; gIndex++)
        {
            GROUPPARMS *group = resiParms->groupList[gIndex];
            for (int aIndex = 0; aIndex < group->nAtoms; aIndex++)
            {
                if (strcmp(atomName, group->atomList[aIndex]->atomName) == 0)
                {
                    parms->speciesInfo[i].name = strdup(spname);
                    parms->speciesInfo[i].rIndex = rIndex;
                    parms->speciesInfo[i].gIndex = gIndex;
                    parms->speciesInfo[i].aIndex = aIndex;
                    parms->speciesInfo[i].offset = offset;
                    parms->speciesInfo[i].isResidueOwner = isResidueOwner;
                }
                offset++;
            }
        }
    }
    parms->isResidueOwnerSpecies = isResidueOwnerArr;

    parms->charmmParms->charmmWeights = NULL;

    // Set up function pointers for connectivity energy terms
    parms->bond_fcn = (double (*)(void *, void *, int, int, void *, void *, void *, void *))resBondSorted;
    parms->angle_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resAngleSorted;
    parms->cosangle_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resAngleCosineSorted;
    parms->rebangle_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resAngleRestrainSorted;
    parms->ureybradley_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))NULL;
    parms->torsion_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resTorsionSorted;
    parms->improper_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resImproperSorted;
    parms->cmap_fcn = NULL;
    parms->bpair_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *, double *, double *, void *))resCGBpairSorted;

    if ((parms->excludePotentialTerm & bondMask) != 0) parms->bond_fcn = NULL;
    if ((parms->excludePotentialTerm & angleMask) != 0) parms->angle_fcn = NULL;
    if ((parms->excludePotentialTerm & cosangleMask) != 0) parms->cosangle_fcn = NULL;
    if ((parms->excludePotentialTerm & rebangleMask) != 0) parms->rebangle_fcn = NULL;
    if ((parms->excludePotentialTerm & ureybradleyMask) != 0) parms->ureybradley_fcn = NULL;
    if ((parms->excludePotentialTerm & torsionMask) != 0) parms->torsion_fcn = NULL;
    if ((parms->excludePotentialTerm & improperMask) != 0) parms->improper_fcn = NULL;
    if ((parms->excludePotentialTerm & cMapMask) != 0) parms->cmap_fcn = NULL;
    if ((parms->excludePotentialTerm & nonBondMask) != 0) parms->bpair_fcn = NULL;

    //check if we want to grab gpu parms
    ACCELERATOR *accelerator = NULL;
    accelerator = accelerator_getAccelerator(accelerator);
    if ((accelerator != NULL) && accelerator->itype == GPU_CUDA)
    {
        if (getRank(0) == 0) printf("using gpu martini parms\n");
        GPUCODE(martiniNonBondGPUParms(parms);)
        GPUCODE(martiniBondGPUParms(parms);)
        potential->use_gpu_list = 1;
        potential->neighborTableType = NEIGHBORTABLE_GPU;
        GPUCODE(potential->eval_potential = (void (*)(void *, void *, void *))martiniGPU1;)
    }
    else
    {
        printf("using cpu martini parms\n");
    }

    timestamp("END   Bio Martini Setup");
    return parms;
}

static BIOENERGIES bioEnergiesZero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void martini(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    static int firstTime = 1;
    //STATE *state=sys->collection->state;
    profile(CHARMM_T, START);
    parms->bioEnergies = bioEnergiesZero;
    if (sys->moleculeClass != NULL)
    {
        if (sys->neighbor->lastUpdate == sys->loop) reOrgPairs(sys, parms);

        if ((parms->excludePotentialTerm & nonBondMask) == 0)
        {
            martiniNonBond(sys, parms, e); // Nonbond interaction including LJ and ELE from Reaction Field/
            martiniIntraMoleReaction(sys, parms, e);
        }
        parms->bpair_fcn = NULL;
    }
    else
    {
        if ((parms->excludePotentialTerm & nonBondMask) == 0)
        {
            martiniNonBond(sys, parms, e); // Nonbond interaction including LJ and ELE from Reaction Field/
        }
    }
    charmmConvalent(sys, parms, e); //Calculate all bonded terms
    if (firstTime)
    {
        BIOENERGIES bioEnergies = parms->bioEnergies;
        MPI_Allreduce(&bioEnergies, &parms->bioEnergies, sizeof (BIOENERGIES) / sizeof (double), MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
        if (getRank(0) == 0) printBioEnergies(parms->bioEnergies);
        firstTime = 0;
    }
    profile(CHARMM_T, END);
}

void reOrgPairs(SYSTEM*sys, CHARMMPOT_PARMS *parms)
{
    assert(sys->neighbor->nc == 1); // assume only a single cutoff for pruned neighbor tables. 
    //int IndexRmax = -1;
    //for (int i = 0; i < sys->neighbor->nc; i++) {
    //   if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
    //}
    //if (IndexRmax == -1) 
    //   error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));

    CHARMM_PARMS *charmmParms = parms->charmmParms;
    NPARTICLE *particles = sys->neighbor->particles;
    //double rcut = parms->rmax;

    STATE *state = sys->collection->state;
    gid_type *label = state->label;
    SPECIES **species = state->species;
    double minspan;
    box_get(NULL, MINSPAN, &minspan);
    //double R2cut = 0.25*minspan*minspan;
    BPAIR_CONN **bpairStart[sys->nspecies];
    int bpairListSize[sys->nspecies];
    int atomListSize[sys->nspecies];
    MOLECULECLASS *moleculeClass = sys->moleculeClass;
    for (int i = 0; i < moleculeClass->nMoleculeTypes; i++)
    {
        char *name = moleculeClass->moleculeTypes[i]->ownershipSpecies->name;
        RESI_CONN *resiConn = findResiConnNew(charmmParms, name);
        bpairStart[i] = resiConn->bpairList;
        bpairListSize[i] = resiConn->bpairListSize;
        atomListSize[i] = resiConn->atomListSize;
    }

    int cnt = 0;
    for (unsigned i = 0; i < sys->nlocal; i++)
    {
        unsigned atmI = (unsigned) (label[i] & atmgrpMask);
        gid_type moleI = label[i] & molMask;
        int iTypeIndex = moleculeClass->speciesIndexToMoleculeIndex[species[i]->index];
        MOLECULETYPE *moleculeType = moleculeClass->moleculeTypes[iTypeIndex];

        PAIRS *f = NULL;
        PAIRS *q = NULL;
        PAIRS *pf = NULL;
        PAIRS *pq = NULL;
        PAIRS *p = particles[i].ifirst[0];
        while (p != NULL)
        {
            int j = p->j;
            //int jTypeIndex = moleculeClass->speciesIndexToMoleculeIndex[species[j]->index];  
            int prune = 0;
            if (moleI == (label[j] & molMask))
            {
                if (moleculeType->nSpecies > 1)
                {
                    assert(moleculeType->nMembers > 0);
                    unsigned atmJ = (unsigned) (label[j] & atmgrpMask);
                    for (int k = 0; k < bpairListSize[iTypeIndex]; k++)
                    {
                        BPAIR_CONN **bpairList = bpairStart[iTypeIndex];
                        unsigned excludeI = bpairList[k]->atmI;
                        unsigned excludeJ = bpairList[k]->atmJ;
                        {
                            if (((atmI == excludeI) && (atmJ == excludeJ)) || ((atmJ == excludeI) && (atmI == excludeJ)))
                            {
                                prune = 1;
                                break;
                            }
                        }
                    }
                }
                else prune = 1;
            }
            if (prune == 0)
            {
                if (f == NULL) f = p;
                if ((q != NULL) && (q->ilink != p)) q->ilink = p;
                q = p;
            }
            else
            {
                if (pf == NULL) pf = p;
                if ((pq != NULL) && (pq->ilink != p)) pq->ilink = p;
                pq = p;
            }
            p = p->ilink;
            cnt++;
        }
        if (q != NULL) q->ilink = NULL;
        if (pq != NULL) pq->ilink = NULL;
        particles[i].ifirst[0] = f;
        particles[i].ifirst[1] = pf;
    }
}
