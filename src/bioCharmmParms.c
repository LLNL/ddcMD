#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "utilities.h"
#include "bioCharmmParms.h"
#include "bioCharmmPar.h"
#include "bioCharmmTop.h"
#include "bioCharmmCovalent.h"
#include "ddcMalloc.h"

typedef struct conn_size_str
{
    int resiConnSize;
    int grpSize;
    int atmSize;
    int bondSize;
    int angleSize;
    int torsSize;
    int imprSize;
    int cmapSize;
    //how many to exclude
    int weightedSize;
} CONN_SIZE;

typedef struct four_int_str
{
    int index[4];
} FOUR_INT;

void getConnSize(CHARMM_PARMS* charmmParms, CONN_SIZE* connSizePtr)
{
    CHARMM_TOP* charmmTop = charmmParms->charmmTop;

    int resCount = 0;
    int grpCount = 0;
    int atmCount = 0;
    int bndCount = 0;
    int impCount = 0;
    int weightCount = 0;
    for (int i = 0; i < charmmTop->resiParmSize; ++i)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        int resType = resiParm->resType;
        if (resType == 1)
        {
            resCount += 3;
            grpCount += resiParm->groupListSize * 3;
            atmCount += resiParm->atomListSize * 3 + 3; //NTER has two more atoms and CTER has one more
            bndCount += resiParm->bondListSize * 3 + 2; //NTER has two more bonds and CTER has same bonds.
            impCount += resiParm->imprListSize * 3 + 1; //NTER has no impr and CTER has one impr
            weightCount += (resiParm->species1Size * resiParm->species2Size);
        }
        else
        {
            resCount += 1;
            grpCount += resiParm->groupListSize;
            atmCount += resiParm->atomListSize;
            bndCount += resiParm->bondListSize;
            impCount += resiParm->imprListSize;
            weightCount += (resiParm->species1Size * resiParm->species2Size);
        }
    }

    connSizePtr->resiConnSize = resCount;
    connSizePtr->grpSize = grpCount;
    connSizePtr->atmSize = atmCount;
    connSizePtr->bondSize = bndCount;
    connSizePtr->imprSize = impCount;
    connSizePtr->weightedSize = weightCount;
}

int atmExistList(char* atmName, ONE_ATOM** atmList, int atmListSize)
{
    for (int i = 0; i < atmListSize; ++i)
    {
        if (strcmp(atmName, atmList[i]->atomI) == 0)
        {
            return 1;
        }
    }
    return 0;
}

int atmIDinList(char* atmName, TATOM_PARMS** atmList, int atmListSize)
{
    for (int i = 0; i < atmListSize; ++i)
    {
        if (strcmp(atmName, atmList[i]->atmName) == 0)
        {
            return i;
        }
    }
    if (strcmp(atmName, "+N") == 0) return (atmListSize);
    if (strcmp(atmName, "-C") == 0) return -2;

    printf("Cannot find %s in the atom list.\n", atmName);
    exit(EXIT_FAILURE);
}

MASS_PARMS* atmTypeinList(char* atmName, RESI_CONN* resiConn)
{
    /*
        char* realName;
        if(atmName[0]=='+' || atmName[0]=='-'){
            realName=atmName+1;
        }else{
            realName=atmName;
        }
     */
    if (strcmp(atmName, "+N") == 0)
    {
        return resiConn->nh1MassParm;
    }

    if (strcmp(atmName, "-C") == 0)
    {
        return resiConn->cMassParm;
    }
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        if (strcmp(atmName, resiConn->atomList[i]->atmName) == 0)
        {
            return resiConn->atomList[i]->atmTypePtr;
        }
    }

    printf("Cannot find %s in the atom list for MASS_PARMS.\n", atmName);
    exit(EXIT_FAILURE);
}

MASS_PARMS* atmTypeIDinList(char* atmName, RESI_CONN* resiConn, int* atmIDptr)
{

    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        if (strcmp(atmName, resiConn->atomList[i]->atmName) == 0)
        {
            *atmIDptr = i;
            return resiConn->atomList[i]->atmTypePtr;
        }
    }

    if (strcmp(atmName, "+N") == 0)
    {
        *atmIDptr = resiConn->atomListSize;
        return resiConn->nh1MassParm;
    }

    if (strcmp(atmName, "-C") == 0)
    {
        *atmIDptr = -2;
        return resiConn->cMassParm;
    }

    printf("Cannot find %s in the atom list for MASS_PARMS.\n", atmName);
    exit(EXIT_FAILURE);

}

int nTerAtmTypeChanged(RESI_PARMS* resiParm, char* atmName)
{


    if (strcmp(atmName, "N") == 0)
    {
        return 1;
    }
    if (strcmp(atmName, "HA") == 0)
    {
        return 1;
    }

    if (strcmp(atmName, "HA1") == 0 && resNameDiffer(resiParm->resName, "GLY") == 0)
    {
        return 1;
    }
    if (strcmp(atmName, "HA2") == 0 && resNameDiffer(resiParm->resName, "GLY") == 0)
    {
        return 1;
    }

    //    if(strcmp(atmName, "HA")==0 && resNameDiffer(resiParm->resName, "PRO")==0){
    //        return 1;
    //    }     

    if (strcmp(atmName, "HD1") == 0 && resNameDiffer(resiParm->resName, "PRO") == 0)
    {
        return 1;
    }

    if (strcmp(atmName, "HD2") == 0 && resNameDiffer(resiParm->resName, "PRO") == 0)
    {
        return 1;
    }

    return 0;
}

BOND_PARMS* getBondParms(CHARMM_PAR* charmmPar, MASS_PARMS* atmItype, MASS_PARMS* atmJtype)
{
    char* atmI = atmItype->atmType;
    char* atmJ = atmJtype->atmType;
    for (int i = 0; i < charmmPar->bondParmSize; ++i)
    {
        BOND_PARMS* bondParm = charmmPar->bondParms[i];
        if (strcmp(atmI, bondParm->atmTypeI) == 0 && strcmp(atmJ, bondParm->atmTypeJ) == 0)
        {
            return bondParm;
        }
        if (strcmp(atmJ, bondParm->atmTypeI) == 0 && strcmp(atmI, bondParm->atmTypeJ) == 0)
        {
            return bondParm;
        }
    }

    printf("Cannot find bond parameter between %s and %s.\n", atmI, atmJ);
    exit(EXIT_FAILURE);
}

IMPROPER_PARMS* getImprParms(CHARMM_PAR* charmmPar, MASS_PARMS* atmItype, MASS_PARMS* atmJtype,
                             MASS_PARMS* atmKtype, MASS_PARMS* atmLtype)
{

    char* atmI = atmItype->atmType;
    char* atmJ = atmJtype->atmType;
    char* atmK = atmKtype->atmType;
    char* atmL = atmLtype->atmType;

    for (int i = 0; i < charmmPar->imprParmsize; ++i)
    {
        IMPROPER_PARMS* imprParm = charmmPar->imprParms[i];
        if (strcmp(imprParm->atmTypeJ, "X") == 0 || strcmp(imprParm->atmTypeK, "X") == 0)
        {
            if (strcmp(atmI, imprParm->atmTypeI) == 0 && strcmp(atmL, imprParm->atmTypeL) == 0)
            {
                return imprParm;
            }
            if (strcmp(atmL, imprParm->atmTypeI) == 0 && strcmp(atmI, imprParm->atmTypeL) == 0)
            {
                return imprParm;
            }
        }
        else
        {
            if (strcmp(atmI, imprParm->atmTypeI) == 0 && strcmp(atmJ, imprParm->atmTypeJ) == 0
                    && strcmp(atmK, imprParm->atmTypeK) == 0 && strcmp(atmL, imprParm->atmTypeL) == 0)
            {
                return imprParm;
            }
            if (strcmp(atmL, imprParm->atmTypeI) == 0 && strcmp(atmK, imprParm->atmTypeJ) == 0
                    && strcmp(atmJ, imprParm->atmTypeK) == 0 && strcmp(atmI, imprParm->atmTypeL) == 0)
            {
                return imprParm;
            }
        }
    }

    printf("Cannot find improper parameter %s - %s - %s - %s.\n", atmI, atmJ, atmK, atmL);
    exit(EXIT_FAILURE);

}

int genBond(CHARMM_PARMS* charmmParms)
{
    CONN_SIZE connSize;
    getConnSize(charmmParms, &connSize);

    charmmParms->resiConnSize = connSize.resiConnSize;
    RESI_CONN *resConnStack = ddcMalloc(connSize.resiConnSize * sizeof (RESI_CONN));
    charmmParms->resiConnList = ddcMalloc(connSize.resiConnSize * sizeof (RESI_CONN*));
    for (int i = 0; i < connSize.resiConnSize; ++i)
    {
        charmmParms->resiConnList[i] = &resConnStack[i];
    }

    TGROUP_PARMS** grpPtrs = ddcMalloc(connSize.grpSize * sizeof (TGROUP_PARMS*));
    TATOM_PARMS** atmPtrs = ddcMalloc(connSize.atmSize * sizeof (TATOM_PARMS*));

    BOND_CONN *bondHeap = ddcMalloc(connSize.bondSize * sizeof (BOND_CONN));
    BOND_CONN **bondPtrs = ddcMalloc(connSize.bondSize * sizeof (BOND_CONN*));
    for (int i = 0; i < connSize.bondSize; ++i)
    {
        bondPtrs[i] = &bondHeap[i];
    }

    IMPR_CONN *imprHeap = ddcMalloc(connSize.imprSize * sizeof (IMPR_CONN));
    IMPR_CONN **imprPtrs = ddcMalloc(connSize.imprSize * sizeof (IMPR_CONN*));
    for (int i = 0; i < connSize.imprSize; ++i)
    {
        imprPtrs[i] = &imprHeap[i];
    }

    // ZEROWEIGHTS *weightHeap=ddcMalloc(connSize.weightedSize*sizeof(ZEROWEIGHTS));
    // ZEROWEIGHTS **weightPtrs=ddcMalloc(connSize.weightedSize*sizeof(ZEROWEIGHTS*));
    // for(int i=0; i<connSize.weightedSize; ++i)
    //      weightPtrs[i]=&weightHeap[i];  

    CHARMM_PAR* charmmPar = charmmParms->charmmPar;

    //For std AA
    CHARMM_TOP* charmmTop = charmmParms->charmmTop;
    MASS_PARMS* nh1MassParm = 0;
    MASS_PARMS* cMassParm = 0;

    for (int i = 0; i < charmmTop->massParmSize; ++i)
    {
        MASS_PARMS* massParm = charmmTop->massParms[i];
        if (strcmp(massParm->atmType, "NH1") == 0)
        {
            nh1MassParm = massParm;
        }
        if (strcmp(massParm->atmType, "C") == 0)
        {
            cMassParm = massParm;
        }
    }

    RESI_PARMS* nTerParm = 0;
    RESI_PARMS* cTerParm = 0;
    RESI_PARMS* glynParm = 0;
    RESI_PARMS* pronParm = 0;
    RESI_PARMS* newNTerParm = 0;

    for (int i = 0; i < charmmTop->presParmSize; ++i)
    {
        RESI_PARMS* resiParm = charmmTop->presParms[i];
        if (resNameDiffer(resiParm->resName, "NTER") == 0)
        {
            nTerParm = resiParm;
            /*
                        for(unsigned j=0; j<nTerParm->groupListSize; j++){
                           printf("group id: %d", nTerParm->groupList[j]->grpID);
                           for(unsigned k=0; k<nTerParm->groupList[j]->grpAtomSize; k++){
                              printf("atom name: %s", nTerParm->groupList[j]->grpAtoms[k]->atmName);
                           }
                        }
             */
        }
        if (resNameDiffer(resiParm->resName, "CTER") == 0)
        {
            cTerParm = resiParm;
        }
        if (resNameDiffer(resiParm->resName, "GLYP") == 0)
        {
            glynParm = resiParm;
        }
        if (resNameDiffer(resiParm->resName, "PROP") == 0)
        {
            pronParm = resiParm;
        }

    }

    int count = 0;
    int grpTotalSize = 0;
    int atmTotalSize = 0;
    int bndTotalSize = 0;
    int impTotalSize = 0;
    //    int weightedTotalSize=0;   

    for (int i = 0; i < charmmTop->resiParmSize; ++i)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        RESI_CONN* resiConn = charmmParms->resiConnList[i];

        resiConn->resID = resiParm->resID;
        strcpy(resiConn->resName, resiParm->resName);
        resiConn->resType = resiParm->resType;
        resiConn->charge = resiParm->charge;
        resiConn->cTer = 0;
        resiConn->nTer = 0;
        resiConn->nh1MassParm = nh1MassParm;
        resiConn->cMassParm = cMassParm;

        resiConn->groupListSize = resiParm->groupListSize;
        resiConn->groupList = &grpPtrs[grpTotalSize];
        for (int j = 0; j < resiParm->groupListSize; ++j)
        {
            grpPtrs[grpTotalSize] = resiParm->groupList[j];
            ++grpTotalSize;
        }

        resiConn->atomListSize = resiParm->atomListSize;
        resiConn->atomList = &atmPtrs[atmTotalSize];
        for (int j = 0; j < resiParm->atomListSize; ++j)
        {
            atmPtrs[atmTotalSize] = resiParm->atomList[j];
            ++atmTotalSize;
        }

        resiConn->bondListSize = resiParm->bondListSize;
        resiConn->bondList = &bondPtrs[bndTotalSize];
        for (int j = 0; j < resiParm->bondListSize; ++j)
        {
            int atomI = -10;
            int atomJ = -10;
            MASS_PARMS* atmItype = atmTypeIDinList(resiParm->bondList[j]->atomI, resiConn, &atomI);
            MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->bondList[j]->atomJ, resiConn, &atomJ);
            bondPtrs[bndTotalSize]->atmI = atomI;
            bondPtrs[bndTotalSize]->atmJ = atomJ;
            bondPtrs[bndTotalSize]->bondPtr = getBondParms(charmmPar, atmItype, atmJtype);
            ++bndTotalSize;
        }

        resiConn->imprListSize = resiParm->imprListSize;
        resiConn->imprList = &imprPtrs[impTotalSize];
        for (int j = 0; j < resiParm->imprListSize; ++j)
        {
            int atomI = -10;
            int atomJ = -10;
            int atomK = -10;
            int atomL = -10;
            MASS_PARMS* atmItype = atmTypeIDinList(resiParm->imprList[j]->atomI, resiConn, &atomI);
            MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->imprList[j]->atomJ, resiConn, &atomJ);
            MASS_PARMS* atmKtype = atmTypeIDinList(resiParm->imprList[j]->atomK, resiConn, &atomK);
            MASS_PARMS* atmLtype = atmTypeIDinList(resiParm->imprList[j]->atomL, resiConn, &atomL);

            imprPtrs[impTotalSize]->atmI = atomI;
            imprPtrs[impTotalSize]->atmJ = atomJ;
            imprPtrs[impTotalSize]->atmK = atomK;
            imprPtrs[impTotalSize]->atmL = atomL;

            imprPtrs[impTotalSize]->imprPtr = getImprParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);
            ++impTotalSize;
        }
        //Weights
        /*        resiConn->weightListSize=resiParm->species1Size*resiParm->species2Size;
                resiConn->weightList=&weightPtrs[weightedTotalSize];
                for(int j=0; j<resiParm->species1Size; ++j){
                    int atomI=-10;
                    MASS_PARMS* atmItype=atmTypeIDinList(resiParm->species1List[j]->atomI, resiConn, &atomI);
                    for(int k=0; k<resiParm->species2Size; k++){
                      int atomJ=-10;
                      MASS_PARMS* atmJtype=atmTypeIDinList(resiParm->species2List[k]->atomI, resiConn, &atomJ);
                      weightPtrs[weightedTotalSize]->atmI=atomI;
                      weightPtrs[weightedTotalSize]->atmJ=atomJ;
                      ++weightedTotalSize;                    
                    }
                }
         */
        if (resiParm->resType == 1)
        {
            //NTER
            if (resNameDiffer(resiParm->resName, "GLY") == 0)
            {
                newNTerParm = glynParm;
            }
            else if (resNameDiffer(resiParm->resName, "PRO") == 0)
            {
                newNTerParm = pronParm;
            }
            else
            {
                newNTerParm = nTerParm;
            }

            resiConn = charmmParms->resiConnList[charmmTop->resiParmSize + count];
            ++count;

            resiConn->resID = resiParm->resID;
            strcpy(resiConn->resName, resiParm->resName);
            resiConn->resType = resiParm->resType;
            resiConn->charge = resiParm->charge;
            resiConn->nTer = 1;
            resiConn->cTer = 0;
            resiConn->nh1MassParm = nh1MassParm;
            resiConn->cMassParm = cMassParm;

            resiConn->groupListSize = resiParm->groupListSize;
            resiConn->groupList = &grpPtrs[grpTotalSize];
            grpPtrs[grpTotalSize] = newNTerParm->groupList[0];
            ++grpTotalSize;
            for (int j = 1; j < resiParm->groupListSize; ++j)
            {
                grpPtrs[grpTotalSize] = resiParm->groupList[j];
                ++grpTotalSize;
            }

            resiConn->atomList = &atmPtrs[atmTotalSize];
            int atmListSize = 0;
            for (int j = 0; j < resiConn->groupListSize; ++j)
            {
                TGROUP_PARMS* grpParm = resiConn->groupList[j];
                atmListSize = atmListSize + grpParm->grpAtomSize;
                for (int k = 0; k < grpParm->grpAtomSize; ++k)
                {
                    atmPtrs[atmTotalSize] = grpParm->grpAtoms[k];
                    ++atmTotalSize;
                }
            }
            resiConn->atomListSize = atmListSize;

            resiConn->bondList = &bondPtrs[bndTotalSize];
            int bondCount = 0;
            for (int j = 0; j < resiParm->bondListSize; ++j)
            {
                if (atmExistList(resiParm->bondList[j]->atomI, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->bondList[j]->atomJ, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1)
                {
                    continue;
                }
                int atomI = -10;
                int atomJ = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(resiParm->bondList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->bondList[j]->atomJ, resiConn, &atomJ);
                bondPtrs[bndTotalSize]->atmI = atomI;
                bondPtrs[bndTotalSize]->atmJ = atomJ;
                bondPtrs[bndTotalSize]->bondPtr = getBondParms(charmmPar, atmItype, atmJtype);
                ++bndTotalSize;
                ++bondCount;
            }

            for (int j = 0; j < newNTerParm->bondListSize; ++j)
            {
                int atomI = -10;
                int atomJ = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(newNTerParm->bondList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(newNTerParm->bondList[j]->atomJ, resiConn, &atomJ);
                bondPtrs[bndTotalSize]->atmI = atomI;
                bondPtrs[bndTotalSize]->atmJ = atomJ;

                bondPtrs[bndTotalSize]->bondPtr = getBondParms(charmmPar, atmItype, atmJtype);
                ++bndTotalSize;
                ++bondCount;
            }
            resiConn->bondListSize = bondCount; //NTER has two more bonds

            resiConn->imprList = &imprPtrs[impTotalSize];
            int imprCount = 0;
            for (int j = 0; j < resiParm->imprListSize; ++j)
            {
                if (atmExistList(resiParm->imprList[j]->atomI, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomJ, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomK, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomL, newNTerParm->delAtmList, newNTerParm->delAtmListSize) == 1)
                {
                    continue;
                }
                // if PRO NTER residue skip N -C CA CD
                if (resNameDiffer(resiConn->resName, "PRO") == 0 && strcmp(resiParm->imprList[j]->atomI, "N") == 0)
                {
                    continue;
                }
                int atomI = -10;
                int atomJ = -10;
                int atomK = -10;
                int atomL = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(resiParm->imprList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->imprList[j]->atomJ, resiConn, &atomJ);
                MASS_PARMS* atmKtype = atmTypeIDinList(resiParm->imprList[j]->atomK, resiConn, &atomK);
                MASS_PARMS* atmLtype = atmTypeIDinList(resiParm->imprList[j]->atomL, resiConn, &atomL);

                imprPtrs[impTotalSize]->atmI = atomI;
                imprPtrs[impTotalSize]->atmJ = atomJ;
                imprPtrs[impTotalSize]->atmK = atomK;
                imprPtrs[impTotalSize]->atmL = atomL;

                imprPtrs[impTotalSize]->imprPtr = getImprParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);
                ++impTotalSize;
                ++imprCount;
            }
            resiConn->imprListSize = imprCount;

            //CTER
            resiConn = charmmParms->resiConnList[charmmTop->resiParmSize + count];
            ++count;

            resiConn->resID = resiParm->resID;
            strcpy(resiConn->resName, resiParm->resName);
            resiConn->resType = resiParm->resType;
            resiConn->charge = resiParm->charge;
            resiConn->nTer = 0;
            resiConn->cTer = 1;
            resiConn->nh1MassParm = nh1MassParm;
            resiConn->cMassParm = cMassParm;

            resiConn->groupListSize = resiParm->groupListSize;
            resiConn->groupList = &grpPtrs[grpTotalSize];
            for (int j = 0; j < resiParm->groupListSize - 1; ++j)
            {
                grpPtrs[grpTotalSize] = resiParm->groupList[j];
                ++grpTotalSize;
            }
            grpPtrs[grpTotalSize] = cTerParm->groupList[0];
            ++grpTotalSize;

            resiConn->atomList = &atmPtrs[atmTotalSize];
            atmListSize = 0;
            for (int j = 0; j < resiConn->groupListSize; ++j)
            {
                TGROUP_PARMS* grpParm = resiConn->groupList[j];
                atmListSize = atmListSize + grpParm->grpAtomSize;
                for (int k = 0; k < grpParm->grpAtomSize; ++k)
                {
                    atmPtrs[atmTotalSize] = grpParm->grpAtoms[k];
                    ++atmTotalSize;
                }
            }
            resiConn->atomListSize = atmListSize;

            resiConn->bondList = &bondPtrs[bndTotalSize];
            bondCount = 0;
            for (int j = 0; j < resiParm->bondListSize; ++j)
            {
                if (atmExistList(resiParm->bondList[j]->atomI, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->bondList[j]->atomJ, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1
                        || strcmp(resiParm->bondList[j]->atomI, "+N") == 0
                        || strcmp(resiParm->bondList[j]->atomJ, "+N") == 0)
                {
                    continue;
                }

                int atomI = -10;
                int atomJ = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(resiParm->bondList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->bondList[j]->atomJ, resiConn, &atomJ);
                bondPtrs[bndTotalSize]->atmI = atomI;
                bondPtrs[bndTotalSize]->atmJ = atomJ;

                bondPtrs[bndTotalSize]->bondPtr = getBondParms(charmmPar, atmItype, atmJtype);
                ++bndTotalSize;
                ++bondCount;

            }

            for (int j = 0; j < cTerParm->bondListSize; ++j)
            {
                int atomI = -10;
                int atomJ = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(cTerParm->bondList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(cTerParm->bondList[j]->atomJ, resiConn, &atomJ);
                bondPtrs[bndTotalSize]->atmI = atomI;
                bondPtrs[bndTotalSize]->atmJ = atomJ;
                bondPtrs[bndTotalSize]->bondPtr = getBondParms(charmmPar, atmItype, atmJtype);
                ++bndTotalSize;
                ++bondCount;
            }
            resiConn->bondListSize = bondCount; //CTER has same number of bonds

            resiConn->imprList = &imprPtrs[impTotalSize];
            imprCount = 0;
            for (int j = 0; j < resiParm->imprListSize; ++j)
            {
                if (atmExistList(resiParm->imprList[j]->atomI, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomJ, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomK, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1
                        || atmExistList(resiParm->imprList[j]->atomL, cTerParm->delAtmList, cTerParm->delAtmListSize) == 1)
                {
                    continue;
                }

                int atomI = -10;
                int atomJ = -10;
                int atomK = -10;
                int atomL = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(resiParm->imprList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(resiParm->imprList[j]->atomJ, resiConn, &atomJ);
                MASS_PARMS* atmKtype = atmTypeIDinList(resiParm->imprList[j]->atomK, resiConn, &atomK);
                MASS_PARMS* atmLtype = atmTypeIDinList(resiParm->imprList[j]->atomL, resiConn, &atomL);

                imprPtrs[impTotalSize]->atmI = atomI;
                imprPtrs[impTotalSize]->atmJ = atomJ;
                imprPtrs[impTotalSize]->atmK = atomK;
                imprPtrs[impTotalSize]->atmL = atomL;

                imprPtrs[impTotalSize]->imprPtr = getImprParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);
                ++impTotalSize;
                ++imprCount;
            }
            for (int j = 0; j < cTerParm->imprListSize; ++j)
            {

                int atomI = -10;
                int atomJ = -10;
                int atomK = -10;
                int atomL = -10;
                MASS_PARMS* atmItype = atmTypeIDinList(cTerParm->imprList[j]->atomI, resiConn, &atomI);
                MASS_PARMS* atmJtype = atmTypeIDinList(cTerParm->imprList[j]->atomJ, resiConn, &atomJ);
                MASS_PARMS* atmKtype = atmTypeIDinList(cTerParm->imprList[j]->atomK, resiConn, &atomK);
                MASS_PARMS* atmLtype = atmTypeIDinList(cTerParm->imprList[j]->atomL, resiConn, &atomL);

                imprPtrs[impTotalSize]->atmI = atomI;
                imprPtrs[impTotalSize]->atmJ = atomJ;
                imprPtrs[impTotalSize]->atmK = atomK;
                imprPtrs[impTotalSize]->atmL = atomL;

                imprPtrs[impTotalSize]->imprPtr = getImprParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);
                ++impTotalSize;
                ++imprCount;
            }
            resiConn->imprListSize = imprCount;
        }
    }

    return 0;

}

int isBonded(BOND_CONN** bondList, int bondListSize, int atomI, int atomJ)
{
    if (atomI > atomJ)
    {
        int atom = atomI;
        atomI = atomJ;
        atomJ = atom;
    }
    int i = 0;
    for (; i < bondListSize; ++i) if (bondList[i]->atmI == atomI) break;
    for (; i < bondListSize; ++i)
    {
        if (bondList[i]->atmI > atomI) break;
        if (bondList[i]->atmJ == atomJ) return 1;
    }
    return 0;
}

int isBondedNew(RESI_CONN *resiConn, int atomI, int atomJ)
{
    BOND_CONN *bondList = *resiConn->bondList;
    int *firstBond = resiConn->firstBond;
    if (atomI > atomJ)
    {
        int atom = atomI;
        atomI = atomJ;
        atomJ = atom;
    }
    for (int ii = firstBond[atomI]; ii < firstBond[atomI + 1]; ii++)
        if (bondList[ii].atmJ == atomJ) return 1;
    return 0;
}

ANGLE_PARMS* getAngleParms(CHARMM_PAR* charmmPar, MASS_PARMS* atmItype, MASS_PARMS* atmJtype,
                           MASS_PARMS* atmKtype)
{

    char* atmI = atmItype->atmType;
    char* atmJ = atmJtype->atmType;
    char* atmK = atmKtype->atmType;

    for (int i = 0; i < charmmPar->angleParmSize; ++i)
    {
        ANGLE_PARMS* angleParm = charmmPar->angleParms[i];

        if (strcmp(atmI, angleParm->atmTypeI) == 0 && strcmp(atmJ, angleParm->atmTypeJ) == 0
                && strcmp(atmK, angleParm->atmTypeK) == 0)
        {
            return angleParm;
        }
        if (strcmp(atmK, angleParm->atmTypeI) == 0 && strcmp(atmJ, angleParm->atmTypeJ) == 0
                && strcmp(atmI, angleParm->atmTypeK) == 0)
        {
            return angleParm;
        }

    }

    printf("Cannot find angle parameter %s - %s - %s.\n", atmI, atmJ, atmK);
    exit(EXIT_FAILURE);

}

int genAngle(CHARMM_PARMS* charmmParms)
{

    int angleTotalSize = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        for (int i = 0; i < resiConn->atomListSize + 1; ++i)
        {
            for (int j = i + 1; j < resiConn->atomListSize + 1; ++j)
            {
                for (int k = j + 1; k < resiConn->atomListSize + 1; ++k)
                {
                    if (isBonded(resiConn->bondList, resiConn->bondListSize, i, j) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, j, k) == 1)
                    {
                        ++angleTotalSize;
                    }
                    if (isBonded(resiConn->bondList, resiConn->bondListSize, j, i) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, i, k) == 1)
                    {
                        ++angleTotalSize;
                    }
                    if (isBonded(resiConn->bondList, resiConn->bondListSize, i, k) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, k, j) == 1)
                    {
                        ++angleTotalSize;
                    }
                }
            }
        }
    }

    // Need additional two angles for each resiConn with cter==1 and neter==0
    angleTotalSize = angleTotalSize + charmmParms->resiConnSize * 4; // It should be more than enough

    ANGLE_CONN *angleHeap = ddcMalloc(angleTotalSize * sizeof (ANGLE_CONN));
    ANGLE_CONN **anglePtrs = ddcMalloc(angleTotalSize * sizeof (ANGLE_CONN*));
    for (int i = 0; i < angleTotalSize; ++i)
    {
        anglePtrs[i] = &angleHeap[i];
    }

    CHARMM_PAR* charmmPar = charmmParms->charmmPar;
    CHARMM_TOP* charmmTop = charmmParms->charmmTop;
    //For std AA
    MASS_PARMS* nh1MassParm = 0;
    MASS_PARMS* cMassParm = 0;
    MASS_PARMS* hMassParm = 0;
    MASS_PARMS* ct1MassParm = 0;
    //For GLY
    MASS_PARMS* ct2MassParm = 0;
    // For PRO
    MASS_PARMS* nMassParm = 0;
    MASS_PARMS* cp3MassParm = 0;
    MASS_PARMS* cp1MassParm = 0;

    for (int i = 0; i < charmmTop->massParmSize; ++i)
    {
        MASS_PARMS* massParm = charmmTop->massParms[i];
        if (strcmp(massParm->atmType, "NH1") == 0)
        {
            nh1MassParm = massParm;
        }
        if (strcmp(massParm->atmType, "C") == 0)
        {
            cMassParm = massParm;
        }
        if (strcmp(massParm->atmType, "H") == 0)
        {
            hMassParm = massParm;
        }
        if (strcmp(massParm->atmType, "CT1") == 0)
        {
            ct1MassParm = massParm;
        }
        if (strcmp(massParm->atmType, "CT2") == 0)
        {
            ct2MassParm = massParm;
        }
        if (strcmp(massParm->atmType, "N") == 0)
        {
            nMassParm = massParm;
        }
        if (strcmp(massParm->atmType, "CP3") == 0)
        {
            cp3MassParm = massParm;
        }
        if (strcmp(massParm->atmType, "CP1") == 0)
        {
            cp1MassParm = massParm;
        }
    }

    ANGLE_PARMS* c_nh1_h = NULL;
    ANGLE_PARMS* c_nh1_ct1 = NULL;
    ANGLE_PARMS* c_nh1_ct2 = NULL;
    ANGLE_PARMS* c_n_cp1 = NULL;
    ANGLE_PARMS* c_n_cp3 = NULL;

    if (cMassParm != NULL && nh1MassParm != NULL && hMassParm != NULL)
        c_nh1_h = getAngleParms(charmmPar, cMassParm, nh1MassParm, hMassParm);

    if (cMassParm != NULL && nh1MassParm != NULL && ct1MassParm != NULL)
        c_nh1_ct1 = getAngleParms(charmmPar, cMassParm, nh1MassParm, ct1MassParm);

    if (cMassParm != NULL && nh1MassParm != NULL && ct2MassParm != NULL)
        c_nh1_ct2 = getAngleParms(charmmPar, cMassParm, nh1MassParm, ct2MassParm);

    if (cMassParm != NULL && nMassParm != NULL && cp1MassParm != NULL)
        c_n_cp1 = getAngleParms(charmmPar, cMassParm, nMassParm, cp1MassParm);

    if (cMassParm != NULL && nMassParm != NULL && cp3MassParm != NULL)
        c_n_cp3 = getAngleParms(charmmPar, cMassParm, nMassParm, cp3MassParm);

    angleTotalSize = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        int angleSize = 0;
        resiConn->angleList = &anglePtrs[angleTotalSize];
        resiConn->cosangleListSize = 0; // CHARMM doesn't use cosine angle

        // For TIP3 water only HT-OT-HT angle counts.
        if (resNameDiffer(resiConn->resName, "TIP3") == 0)
        {
            if (resiConn->atomListSize != 3)
            {
                printf("Error: TIP3 residue has %d atoms instead of 3\n", resiConn->atomListSize);
                exit(1);
            }
            anglePtrs[angleTotalSize]->atmI = 1; // H1
            anglePtrs[angleTotalSize]->atmJ = 0; // OH2
            anglePtrs[angleTotalSize]->atmK = 2; // H2           
            MASS_PARMS* atmItype = resiConn->atomList[1]->atmTypePtr; // H1
            MASS_PARMS* atmJtype = resiConn->atomList[0]->atmTypePtr; // OH2
            MASS_PARMS* atmKtype = resiConn->atomList[2]->atmTypePtr; // H2
            anglePtrs[angleTotalSize]->anglePtr = getAngleParms(charmmPar, atmItype, atmJtype, atmKtype);
            ++angleTotalSize;
            ++angleSize;
            resiConn->angleListSize = angleSize;
            continue; // end of angle conn calculation for TIP3.
        }

        // within the residue
        for (int i = 0; i < resiConn->atomListSize - 1; ++i)
        {
            for (int j = i + 1; j < resiConn->atomListSize; ++j)
            {
                for (int k = j + 1; k < resiConn->atomListSize + 1; ++k)
                {
                    if (isBonded(resiConn->bondList, resiConn->bondListSize, i, j) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, j, k) == 1)
                    {
                        anglePtrs[angleTotalSize]->atmI = i;
                        anglePtrs[angleTotalSize]->atmJ = j;
                        anglePtrs[angleTotalSize]->atmK = k;

                        MASS_PARMS* atmItype = resiConn->atomList[i]->atmTypePtr;
                        MASS_PARMS* atmJtype = resiConn->atomList[j]->atmTypePtr;
                        MASS_PARMS* atmKtype;
                        if (k == resiConn->atomListSize)
                        {
                            atmKtype = nh1MassParm;
                        }
                        else
                        {
                            atmKtype = resiConn->atomList[k]->atmTypePtr;
                        }

                        //printf("Residue Name %s, Nter=%d, Cter=%d,  %d - %d - %d.\n", 
                        //        resiConn->resName, resiConn->nTer, resiConn->cTer, i, j, k);  //resClean

                        anglePtrs[angleTotalSize]->anglePtr = getAngleParms(charmmPar, atmItype, atmJtype, atmKtype);

                        ++angleTotalSize;
                        ++angleSize;
                    }

                    if (isBonded(resiConn->bondList, resiConn->bondListSize, j, i) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, i, k) == 1)
                    {
                        anglePtrs[angleTotalSize]->atmI = j;
                        anglePtrs[angleTotalSize]->atmJ = i;
                        anglePtrs[angleTotalSize]->atmK = k;

                        MASS_PARMS* atmItype = resiConn->atomList[i]->atmTypePtr;
                        MASS_PARMS* atmJtype = resiConn->atomList[j]->atmTypePtr;
                        MASS_PARMS* atmKtype;
                        if (k == resiConn->atomListSize)
                        {
                            atmKtype = nh1MassParm;
                        }
                        else
                        {
                            atmKtype = resiConn->atomList[k]->atmTypePtr;
                        }

                        //printf("Residue Name %s, Nter=%d, Cter=%d,  %d - %d - %d.\n", 
                        //        resiConn->resName, resiConn->nTer, resiConn->cTer, j, i, k);  //resClean

                        anglePtrs[angleTotalSize]->anglePtr = getAngleParms(charmmPar, atmJtype, atmItype, atmKtype);

                        ++angleTotalSize;
                        ++angleSize;
                    }
                    if (isBonded(resiConn->bondList, resiConn->bondListSize, i, k) == 1
                            && isBonded(resiConn->bondList, resiConn->bondListSize, k, j) == 1)
                    {
                        anglePtrs[angleTotalSize]->atmI = i;
                        anglePtrs[angleTotalSize]->atmJ = k;
                        anglePtrs[angleTotalSize]->atmK = j;

                        MASS_PARMS* atmItype = resiConn->atomList[i]->atmTypePtr;
                        MASS_PARMS* atmJtype = resiConn->atomList[j]->atmTypePtr;
                        MASS_PARMS* atmKtype;
                        if (k == resiConn->atomListSize)
                        {
                            atmKtype = nh1MassParm;
                        }
                        else
                        {
                            atmKtype = resiConn->atomList[k]->atmTypePtr;
                        }
                        //printf("Residue Name %s, Nter=%d, Cter=%d,  %d - %d - %d.\n", 
                        //        resiConn->resName, resiConn->nTer, resiConn->cTer, i, k, j); //resClean

                        anglePtrs[angleTotalSize]->anglePtr = getAngleParms(charmmPar, atmItype, atmKtype, atmJtype);

                        ++angleTotalSize;
                        ++angleSize;
                    }
                }
            }
        }

        //between the residue only for std AA and nter==0
        if (resiConn->resType == 1)
        {
            if (resiConn->nTer == 0)
            {
                // if it is proline treat differently
                if (resNameDiffer(resiConn->resName, "PRO") == 0)
                {
                    anglePtrs[angleTotalSize]->atmI = 1;
                    anglePtrs[angleTotalSize]->atmJ = 0;
                    anglePtrs[angleTotalSize]->atmK = -2; // Put C into atmK avoid crash in genAtmRange
                    anglePtrs[angleTotalSize]->anglePtr = c_n_cp3;
                    ++angleTotalSize;
                    ++angleSize;
                    anglePtrs[angleTotalSize]->atmI = 4;
                    anglePtrs[angleTotalSize]->atmJ = 0;
                    anglePtrs[angleTotalSize]->atmK = -2; // Put C into atmK avoid crash in genAtmRange
                    anglePtrs[angleTotalSize]->anglePtr = c_n_cp1;
                    ++angleTotalSize;
                    ++angleSize;
                }
                else
                {
                    anglePtrs[angleTotalSize]->atmI = 1;
                    anglePtrs[angleTotalSize]->atmJ = 0;
                    anglePtrs[angleTotalSize]->atmK = -2; // Put C into atmK avoid crash in genAtmRange
                    anglePtrs[angleTotalSize]->anglePtr = c_nh1_h;
                    ++angleTotalSize;
                    ++angleSize;
                    anglePtrs[angleTotalSize]->atmI = 2;
                    anglePtrs[angleTotalSize]->atmJ = 0;
                    anglePtrs[angleTotalSize]->atmK = -2; // Put C into atmK avoid crash in genAtmRange 
                    if (resNameDiffer(resiConn->resName, "GLY") == 0)
                    {
                        anglePtrs[angleTotalSize]->anglePtr = c_nh1_ct2;
                    }
                    else
                    {
                        anglePtrs[angleTotalSize]->anglePtr = c_nh1_ct1;
                    }
                    ++angleTotalSize;
                    ++angleSize;
                }
            }
        }

        resiConn->angleListSize = angleSize;
    }

    return 0;
}

TORSION_PARMS* getTorsParms(CHARMM_PAR* charmmPar, MASS_PARMS* atmItype, MASS_PARMS* atmJtype,
                            MASS_PARMS* atmKtype, MASS_PARMS* atmLtype)
{

    char* atmI = atmItype->atmType;
    char* atmJ = atmJtype->atmType;
    char* atmK = atmKtype->atmType;
    char* atmL = atmLtype->atmType;

    for (int i = 0; i < charmmPar->torParmSize; ++i)
    {

        TORSION_PARMS* torsParm = charmmPar->torParms[i];
        if (strcmp("X", torsParm->atmTypeI) == 0 || strcmp("X", torsParm->atmTypeL) == 0)
        {
            if (strcmp(atmJ, torsParm->atmTypeJ) == 0
                    && strcmp(atmK, torsParm->atmTypeK) == 0)
            {
                return torsParm;
            }
            if (strcmp(atmK, torsParm->atmTypeJ) == 0
                    && strcmp(atmJ, torsParm->atmTypeK) == 0)
            {
                return torsParm;
            }
        }
        else
        {
            if (strcmp(atmI, torsParm->atmTypeI) == 0 && strcmp(atmJ, torsParm->atmTypeJ) == 0
                    && strcmp(atmK, torsParm->atmTypeK) == 0 && strcmp(atmL, torsParm->atmTypeL) == 0)
            {
                return torsParm;
            }
            if (strcmp(atmL, torsParm->atmTypeI) == 0 && strcmp(atmK, torsParm->atmTypeJ) == 0
                    && strcmp(atmJ, torsParm->atmTypeK) == 0 && strcmp(atmI, torsParm->atmTypeL) == 0)
            {
                return torsParm;
            }

        }

    }

    printf("Cannot find torsion parameter %s - %s - %s - %s.\n", atmI, atmJ, atmK, atmL);
    exit(EXIT_FAILURE);

}

MASS_PARMS* getMassParmsTors(int atomID, RESI_CONN* resiConn, MASS_PARMS** interMassParm)
{
    if (atomID >= 0)
    {
        if (atomID < resiConn->atomListSize)
        {
            return resiConn->atomList[atomID]->atmTypePtr;
        }
        else if (atomID < resiConn->atomListSize + 3)
        {
            return interMassParm[atomID - resiConn->atomListSize + 1];
        }
        else
        {
            printf("Atom ID out of scope  %d.\n", atomID);
            exit(EXIT_FAILURE);
        }
    }
    // for -C -2 
    return interMassParm[0];

}

int genTorsion(CHARMM_PARMS* charmmParms)
{
    /*  
        FOUR_INT indexPermu[12];
        int count=0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (j != i) {
                    for (int k = 0; k < 4; ++k) {
                        if (k != i && k != j) {
                            for (int l = 0; l < 4; ++l) {
                                if (l != i && l != j && l != k) {
                                    //printf("%d%d%d%d\n", i, j, k, l);
                                    indexPermu[count].index[0]=i;
                                    indexPermu[count].index[1]=j;
                                    indexPermu[count].index[2]=k;
                                    indexPermu[count].index[3]=l;
                                    ++count;
                                }
                            }
                        }
                    }
                }
            }
        } 
     */

    int torsTotalSize = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        //int *first = resiConn->firstBond; 
        /*
            for(int i=0; i<resiConn->atomListSize+1; ++i){
                for(int j=i+1; j<resiConn->atomListSize+1; ++j){
                    for(int k=j+1; k<resiConn->atomListSize+1; ++k){
                        for(int l=k+1; l<resiConn->atomListSize+1; ++l){
                            int atm[4];
                            atm[0]=i;
                            atm[1]=j;
                            atm[2]=k;
                            atm[3]=l;
                            for(int index=0; index<12; ++index){
                                int indexI=indexPermu[index].index[0];
                                int indexJ=indexPermu[index].index[1];
                                int indexK=indexPermu[index].index[2];
                                int indexL=indexPermu[index].index[3];
                                int atmI=atm[indexI];
                                int atmJ=atm[indexJ];
                                int atmK=atm[indexK];
                                int atmL=atm[indexL];
                            
                                if(   isBonded(resiConn->bondList+0*first[atmI], resiConn->bondListSize-0*first[atmI], atmI, atmJ)==1
                                   && isBonded(resiConn->bondList+0*first[atmJ], resiConn->bondListSize-0*first[atmJ], atmJ, atmK)==1 
                                   && isBonded(resiConn->bondList+0*first[atmK], resiConn->bondListSize-0*first[atmK], atmK, atmL)==1 ){
                                    ++torsTotalSize;
                                } 

                            }
                        }
                    }                
                }
            }
         */
        int torsResSize = 0;
        for (int ii = 0; ii < resiConn->bondListSize; ++ii)
        {
            int atmI = resiConn->bondList[ii]->atmI;
            int atmJ = resiConn->bondList[ii]->atmJ;
            torsResSize += (resiConn->bondCnt[atmI] - 1)*(resiConn->bondCnt[atmJ] - 1);
            torsTotalSize += (resiConn->bondCnt[atmI] - 1)*(resiConn->bondCnt[atmJ] - 1); // This will over count torsional terms, because cases where the two end bonds may share an atom. 
        }
        //printf("resname=%s, nter=%d, cter=%d, numTor=%d\n", resiConn->resName, resiConn->nTer, resiConn->cTer, torsResSize);   //resClean
    }

    TORS_CONN *torsHeap = ddcMalloc(torsTotalSize * sizeof (TORS_CONN));
    TORS_CONN **torsPtrs = ddcMalloc(torsTotalSize * sizeof (TORS_CONN*));
    for (int i = 0; i < torsTotalSize; ++i)
    {
        torsPtrs[i] = &torsHeap[i];
    }

    CHARMM_PAR* charmmPar = charmmParms->charmmPar;
    CHARMM_TOP* charmmTop = charmmParms->charmmTop;
    // A help mass array to figure out mass parms for atoms not in current residue.
    MASS_PARMS** interMassParm = ddcMalloc(INTERMASSPNUM * sizeof (MASS_PARMS*));
    for (int i = 0; i < charmmTop->massParmSize; ++i)
    {
        MASS_PARMS* massParm = charmmTop->massParms[i];
        if (strcmp(massParm->atmType, "C") == 0)
        {
            // for -2
            interMassParm[0] = massParm;
        }
        if (strcmp(massParm->atmType, "NH1") == 0)
        {
            // for atomListSize+1
            interMassParm[1] = massParm;
        }
        if (strcmp(massParm->atmType, "H") == 0)
        {
            // for atomListSize+2
            interMassParm[2] = massParm;
        }
        if (strcmp(massParm->atmType, "CT1") == 0)
        {
            // for atomListSize+3
            interMassParm[3] = massParm;
        }
    }

    torsTotalSize = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        int torsSize = 0;
        resiConn->torsList = &torsPtrs[torsTotalSize];

        //TIP3 water residue would have artificial dihedral in Shake. H1-H2 bond.
        if (resNameDiffer(resiConn->resName, "TIP3") == 0)
        {
            resiConn->torsListSize = 0;
            continue;
        }

        for (int b = 0; b < resiConn->bondListSize; ++b)
        {
            int atmI = resiConn->bondList[b]->atmI;
            int atmJ = resiConn->bondList[b]->atmJ;
            MASS_PARMS* atmJtype = getMassParmsTors(atmI, resiConn, interMassParm);
            MASS_PARMS* atmKtype = getMassParmsTors(atmJ, resiConn, interMassParm);
            //torsSize += (resiConn->bondCnt[atmI]-1)*(resiConn->bondCnt[atmJ]-1);
            //torsTotalSize += (resiConn->bondCnt[atmI]-1)*(resiConn->bondCnt[atmJ]-1);   // This will over count torsional terms, because cases where the two end bonds may share an atom. 
            for (int ii = 0; ii < resiConn->bondCnt[atmI]; ii++)
            {
                int torsAtmI = resiConn->bndCntList[atmI].atom[ii];
                if (torsAtmI == atmJ)
                {
                    //If the bond atom is atmJ skip it
                    continue;
                }
                MASS_PARMS* atmItype = getMassParmsTors(torsAtmI, resiConn, interMassParm);
                for (int jj = 0; jj < resiConn->bondCnt[atmJ]; jj++)
                {
                    int torsAtmL = resiConn->bndCntList[atmJ].atom[jj];
                    if (torsAtmL == atmI)
                    {
                        //If the bond atom is atmI skip it
                        continue;
                    }
                    MASS_PARMS* atmLtype = getMassParmsTors(torsAtmL, resiConn, interMassParm);
                    torsPtrs[torsTotalSize]->atmI = torsAtmI;
                    torsPtrs[torsTotalSize]->atmJ = atmI;
                    torsPtrs[torsTotalSize]->atmK = atmJ;
                    torsPtrs[torsTotalSize]->atmL = torsAtmL;

                    torsPtrs[torsTotalSize]->torsPtr = getTorsParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);

                    torsSize++;
                    torsTotalSize++;

                }
            }

        }

        resiConn->torsListSize = torsSize;
    }

    /*
    for(int r=0; r<charmmParms->resiConnSize; ++r){
        RESI_CONN* resiConn=charmmParms->resiConnList[r];
        int torsSize=0;
        resiConn->torsList=&torsPtrs[torsTotalSize];
        for(int i=0; i<resiConn->atomListSize+1; ++i){
            for(int j=i+1; j<resiConn->atomListSize+1; ++j){
                for(int k=j+1; k<resiConn->atomListSize+1; ++k){
                    for(int l=k+1; l<resiConn->atomListSize; ++l){
                        int atm[4];
                        atm[0]=i;
                        atm[1]=j;
                        atm[2]=k;
                        atm[3]=l;
                        for(int index=0; index<12; ++index){
                            int indexI=indexPermu[index].index[0];
                            int indexJ=indexPermu[index].index[1];
                            int indexK=indexPermu[index].index[2];
                            int indexL=indexPermu[index].index[3];
                            int atmI=atm[indexI];
                            int atmJ=atm[indexJ];
                            int atmK=atm[indexK];
                            int atmL=atm[indexL];
                            
                            //if (!isBonded(resiConn->bondList, resiConn->bondListSize, atmI, atmJ)) continue;
                            //if (!isBonded(resiConn->bondList, resiConn->bondListSize, atmJ, atmK)) continue;
                            //if (!isBonded(resiConn->bondList, resiConn->bondListSize, atmK, atmL)) continue; 
                            
                            if(!isBondedNew(resiConn, atmI, atmJ)) continue;
                            if(!isBondedNew(resiConn, atmJ, atmK)) continue; 
                            if(!isBondedNew(resiConn, atmK, atmL)) continue; 
                                
                                torsPtrs[torsTotalSize]->atmI=atmI;
                                torsPtrs[torsTotalSize]->atmJ=atmJ;
                                torsPtrs[torsTotalSize]->atmK=atmK;
                                torsPtrs[torsTotalSize]->atmL=atmL;

                                MASS_PARMS* atmItype;
                                if(atmI==resiConn->atomListSize){
                                    atmItype=nh1MassParm;
                                }else{
                                    atmItype=resiConn->atomList[atmI]->atmTypePtr;
                                }                                
                                MASS_PARMS* atmJtype;
                                if(atmJ==resiConn->atomListSize){
                                    atmJtype=nh1MassParm;
                                }else{
                                    atmJtype=resiConn->atomList[atmJ]->atmTypePtr;
                                } 
                                MASS_PARMS* atmKtype;
                                if(atmK==resiConn->atomListSize){
                                    atmKtype=nh1MassParm;
                                }else{
                                    atmKtype=resiConn->atomList[atmK]->atmTypePtr;
                                }
                                MASS_PARMS* atmLtype;
                                if(atmL==resiConn->atomListSize){
                                    atmLtype=nh1MassParm;
                                }else{
                                    atmLtype=resiConn->atomList[atmL]->atmTypePtr;
                                }
                                

                                //printf("Residue Name %s, Nter=%d, Cter=%d,  %d - %d - %d - %d.\n", 
                                //        resiConn->resName, resiConn->nTer, resiConn->cTer,         //resClean
                                //        atmI, atmJ, atmK, atmL);


                                torsPtrs[torsTotalSize]->torsPtr=getTorsParms(charmmPar, atmItype, atmJtype, atmKtype, atmLtype);

                                //Save atmType for BondPair calculation
                                //strcpy(torsPtrs[torsTotalSize]->atmTypeI, atmItype->atmType);
                                //strcpy(torsPtrs[torsTotalSize]->atmTypeL, atmLtype->atmType);
                                
                                ++torsSize;                                                                
                                ++torsTotalSize;

                        }
                    }
                }                
            }
        }
        resiConn->torsListSize=torsSize;
    }
    
     */

    return 0;
}

#include <assert.h>
void biospline(double dx, double* y, int N, double *u, double *y2) { assert(0);}   //This capability is not support in this release See full version
void splineInterpolation(double xmin, double dx, double *ya, double *y2a, double x, double *y, double *y1){assert(0);}   //This capability is not support in this release. See full version

void computeCMAPDerivatives(CMAP_PARMS* cmapParms)
{
    //cmapParms->map_y1, map_y2, map_y12

    int size = 24 + 12 + 12;
    int nVals = 24;
    int stride = 12;
    double *tempMap = (double *) malloc(size * size * sizeof (double));
    double *u = (double*) malloc(size * sizeof (double));
    double *u2 = (double*) malloc(size * sizeof (double));
    double *t = (double*) malloc(size * size * sizeof (double));
    double dx = 15.0;
    double *ytemp = (double *) malloc(size * sizeof (double));
    double *y1temp = (double *) malloc(size * sizeof (double));
    //offsets for degress
    double degreeOffset = 180.0;
    double minDegree = -360.0;

    double v, v1, v2, v12;

    //copy to temp array
    for (int i = 0; i < size; i++)
    {
        int ii = (i + nVals - stride) % nVals;
        for (int j = 0; j < size; j++)
        {
            int jj = (j + nVals - stride) % nVals;
            tempMap[i * size + j] = cmapParms->grid[ii * nVals + jj];
        }
    }
    for (int i = 0; i < size; i++)
    {
        biospline(dx, &(tempMap[i * size]), size, u, &(t[i * size]));
    }
    for (int i = stride; i < nVals + stride; i++)
    {
        double phi = (i - stride) * dx - degreeOffset;
        for (int j = stride; j < (nVals + stride); j++)
        {
            double psi = (j - stride) * dx - degreeOffset;
            for (int k = 0; k < size; k++)
            {
                splineInterpolation(minDegree, dx, &(tempMap[k * size]), &(t[k * size]), psi, &(ytemp[k]), &(y1temp[k]));
            }
            biospline(dx, ytemp, size, u, u2);
            splineInterpolation(minDegree, dx, ytemp, u2, phi, &v, &v1);
            biospline(dx, y1temp, size, u, u2);
            splineInterpolation(minDegree, dx, y1temp, u2, phi, &v2, &v12);
            cmapParms->map_y1[(i - stride) * nVals + (j - stride)] = v1;
            cmapParms->map_y2[(i - stride) * nVals + (j - stride)] = v2;
            cmapParms->map_y12[(i - stride) * nVals + (j - stride)] = v12;
        }
    }
    free(tempMap);
    free(u);
    free(u2);
    free(t);
    free(ytemp);
    free(y1temp);
}

void precomputeCMAP(CHARMM_PARMS* charmmParms)
{

    CHARMM_PAR* charmmPar = charmmParms->charmmPar;
    for (int i = 0; i < charmmPar->cmapParmSize; i++)
    {
        CMAP_PARMS* cmapParms = charmmPar->cmapParms[i];
        computeCMAPDerivatives(cmapParms);
    }
}

int checkForCmapTerm(RESI_CONN* resiConn, CHARMM_TOP* charmmTop)
{
    for (int i = 0; i < charmmTop->resiParmSize; i++)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        if (resNameDiffer(resiConn->resName, resiParm->resName) == 0)
        {
            if (resiParm->cmapListSize == 0) return -1;
            return 0;
        }
    }
    return -1;
}

int getCMAPParm(RESI_CONN* resiConn, CHARMM_TOP* charmmTop, int type, char* atom)
{
    for (int i = 0; i < charmmTop->resiParmSize; i++)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        if (resNameDiffer(resiConn->resName, resiParm->resName) == 0)
        {
            if (resiParm->cmapListSize == 0) return -1;
            if (type == 0)
            {
                strncpy(atom, resiParm->cmapList[1]->atomK, 5);
                return 0;
            }
            if (type == 1)
            {
                strncpy(atom, resiParm->cmapList[0]->atomJ, 5);
                return 0;
            }
            if (type == 2)
            {
                strncpy(atom, resiParm->cmapList[0]->atomK, 5);
                return 0;
            }
            if (type == 3)
            {
                strncpy(atom, resiParm->cmapList[0]->atomL, 5);
                return 0;
            }
            if (type == 4)
            {
                strncpy(atom, resiParm->cmapList[0]->atomJ, 5);
                return 0;
            }
        }
    }
    return -1;
}

int getCMAPParms(RESI_CONN* resiConn, CHARMM_TOP* charmmTop, char* atomI, char* atomJ, char* atomK)
{

    for (int i = 0; i < charmmTop->resiParmSize; i++)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        if (resNameDiffer(resiConn->resName, resiParm->resName) == 0)
        {
            if (resiParm->cmapListSize == 0) return -1;
            strncpy(atomI, resiParm->cmapList[0]->atomJ, 5);
            strncpy(atomJ, resiParm->cmapList[0]->atomK, 5);
            strncpy(atomK, resiParm->cmapList[0]->atomL, 5);
            return 0;
        }
    }
    return -1;
}

int getAtomIds(RESI_CONN* resiConn, char* atom)
{
    //TATOM_PARMS** atomParms=resiConn->atomList;
    for (int i = 0; i < resiConn->atomListSize; i++)
    {
        TATOM_PARMS* atomParms = resiConn->atomList[i];
        if (strncmp(atom, atomParms->atmName, 5) == 0) return i;
    }
    return -1;
}

//Types:
// ALA: 0
// ALA before PRO: 1
// PRO: 2
// 2 Adjacent Prolines: 3
// GLY: 4
// GLY before Pro: 5
// Connectivity: (A)-B-C

int getCmapType(RESI_CONN* resiConn, int prolineFlag)
{
    char resB[6];
    strncpy(resB, resiConn->resName, 6);
    if (strncmp(resB, "ALA", 3) == 0)
    {
        if (prolineFlag) return 1;
        return 0;
    }
    if (strncmp(resB, "PRO", 3) == 0)
    {
        if (prolineFlag) return 3;
        return 2;
    }
    if (strncmp(resB, "GLY", 3) == 0)
    {
        if (prolineFlag) return 5;
        return 4;
    }

    return -1;
}

int genCmap(CHARMM_PARMS* charmmParms)
{

    CHARMM_TOP* charmmTop = charmmParms->charmmTop;

    long int cmapTotalSize = 0;
    //Loop through everything once to count    
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[i];
        char res[6];
        strncpy(res, resiConn->resName, 6);

        if (((strncmp(res, "ALA", 3) == 0 && strncmp(res, "ALAD", 4) != 0) || strncmp(res, "GLY", 3) == 0 || strncmp(res, "PRO", 3) == 0) && !resiConn->nTer && !resiConn->cTer)
        {

            cmapTotalSize += 2;
        }
    }

    CMAP_CONN *cmapHeap = ddcMalloc(cmapTotalSize * sizeof (CMAP_CONN));
    CMAP_CONN **cmapPtrs = ddcMalloc(cmapTotalSize * sizeof (CMAP_CONN*));
    for (int i = 0; i < cmapTotalSize; i++)
    {
        cmapPtrs[i] = &cmapHeap[i];
    }
    //Loop through again to store cmap atom info
    cmapTotalSize = 0;
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[i];
        char res[6];
        strncpy(res, resiConn->resName, 6);
        resiConn->cmapListSize = 0;
        if ((strncmp(res, "ALA", 3) == 0 || strncmp(res, "GLY", 3) == 0 || strncmp(res, "PRO", 3) == 0) && !resiConn->nTer && !resiConn->cTer && !(strncmp(res, "ALAD", 4) == 0))
        {
            int cmapSize = 0;

            resiConn->cmapList = &cmapPtrs[cmapTotalSize];
            //store central residue atoms
            char atomI[5];
            char atomJ[5];
            char atomK[5];
            if (getCMAPParms(resiConn, charmmTop, atomI, atomJ, atomK) < 0)
            {
                printf("Cannot find CMAP parameters for residue %s\n", resiConn->resName);
            }

            int atomi = getAtomIds(resiConn, atomI);
            int atomj = getAtomIds(resiConn, atomJ);
            int atomk = getAtomIds(resiConn, atomK);
            //NOTE: this stores the atom index within a group
            cmapPtrs[cmapTotalSize]->atmI = atomi;
            cmapPtrs[cmapTotalSize]->atmJ = atomj;
            cmapPtrs[cmapTotalSize]->atmK = atomk;

            //Determine cmapType
            int prolineFlag = 0;
            cmapPtrs[cmapTotalSize]->cmapType = getCmapType(resiConn, prolineFlag);
            cmapTotalSize++;
            cmapSize++;

            prolineFlag = 1;
            cmapPtrs[cmapTotalSize]->atmI = atomi;
            cmapPtrs[cmapTotalSize]->atmJ = atomj;
            cmapPtrs[cmapTotalSize]->atmK = atomk;
            cmapPtrs[cmapTotalSize]->cmapType = getCmapType(resiConn, prolineFlag);
            cmapTotalSize++;
            cmapSize++;
            resiConn->cmapListSize = cmapSize;
        }
    } // end second section through resiConn 
    return 0;
}

int genWeights(CHARMM_PARMS* charmmParms)
{

    CONN_SIZE connSize;
    getConnSize(charmmParms, &connSize);
    CHARMM_TOP* charmmTop = charmmParms->charmmTop;

    ZEROWEIGHTS *weightHeap = ddcMalloc(connSize.weightedSize * sizeof (ZEROWEIGHTS));
    ZEROWEIGHTS **weightPtrs = ddcMalloc(connSize.weightedSize * sizeof (ZEROWEIGHTS*));
    for (int i = 0; i < connSize.weightedSize; ++i)
        weightPtrs[i] = &weightHeap[i];

    int weightedTotalSize = 0;
    for (int i = 0; i < charmmTop->resiParmSize; ++i)
    {
        RESI_PARMS* resiParm = charmmTop->resiParms[i];
        RESI_CONN* resiConn = charmmParms->resiConnList[i];

        //Weights
        resiConn->weightListSize = resiParm->species1Size * resiParm->species2Size;
        if (resiConn->weightListSize == 0)
        {
            resiConn->weightList = NULL;
            continue;
        }
        resiConn->weightList = &weightPtrs[weightedTotalSize];
        //TODO: check atom numbering scheme for residues with more than one group
        for (int j = 0; j < resiParm->species1Size; ++j)
        {
            int atomI = -10;
            atmTypeIDinList(resiParm->species1List[j]->atomI, resiConn, &atomI);
            for (int k = 0; k < resiParm->species2Size; k++)
            {
                int atomJ = -10;
                atmTypeIDinList(resiParm->species2List[k]->atomI, resiConn, &atomJ);
                weightPtrs[weightedTotalSize]->atmI = atomI;
                weightPtrs[weightedTotalSize]->atmJ = atomJ;
                ++weightedTotalSize;
            }
        }
    }
    return 0;
}

/*
int genWeights(CHARMM_PARMS* charmmParms)
{
   //count exclusions
    int weightTotalSize=0;
    for(int r=0; r<charmmParms->resiConnSize; ++r){
        RESI_CONN* resiConn=charmmParms->resiConnList[r];
        for(int i=0; i<resiConn->atomListSize+1; ++i){
            TATOM_PARMS* atmI = resiConn->atomList[i];
            for(int j=i+1; j<resiConn->atomListSize+1; ++j){
              TATOM_PARMS* atmJ = resiConn->atomList[i];
              double w = get2Weights(charmmParms->charmmWeights, atmI, atmJ);
              if(w==0.0)
                weightTotalSize++;
            }
        }
    }
    //allocate lists

    WEIGHT_CONN *weightHeap=ddcMalloc(weightTotalSize*sizeof(WEIGHT_CONN));
    WEIGHT_CONN **weightPtrs=ddcMalloc(weightTotalSize*sizeof(WEIGHT_CONN*));
    //generate list
    //Note: this assumes one morphing molecule is allowed     
    for(int r=0; r<charmmParms->resiConnSize; ++r){
        RESI_CONN* resiConn=charmmParms->resiConnList[r];
        int mySize = 0;
        resiConn->weightList=&weightPtrs[weightTotalSize];
        for(int i=0; i<resiConn->atomListSize+1; ++i){
            TATOM_PARMS* atmI = resiConn->atomList[i];
            for(int j=i+1; j<resiConn->atomListSize+1; ++j){
              TATOM_PARMS* atmJ = resiConn->atomList[i];
              double w = get2Weights(charmmParms->charmmWeights, atmI, atmJ);
              if(w==0.0){

                weightPtrs[weightTotalSize]->atmI = i;
                weightPtrs[weightTotalSize]->atmJ = j;
                weightTotalSize++;
                mySize++;
              }
            }
        }
        resiConn->weightListSize = mySize;
    }
}
 */
int getLJParm(CHARMM_PAR* charmmPar, char* atmType, double* eps, double* sigma)
{
    for (int i = 0; i < charmmPar->ljParmSize; ++i)
    {
        LJCH_PARMS* ljParms = charmmPar->ljParms[i];
        if (strcmp(ljParms->atmTypeI, atmType) == 0)
        {
            *eps = ljParms->epsilon;
            *sigma = ljParms->rmin;
            return 0;
        }
    }

    return -1;
}

int getLJpairParm(CHARMM_PAR* charmmPar, char* atmTypeI, char* atmTypeJ, double* eps, double* sigma)
{
    double epsI = 0;
    double epsJ = 0;
    double sigmaI = 0;
    double sigmaJ = 0;
    if (getLJParm(charmmPar, atmTypeI, &epsI, &sigmaI) < 0)
    {
        printf("BondedPair cannot find %s  LJ parameter\n", atmTypeI);
    }
    if (getLJParm(charmmPar, atmTypeJ, &epsJ, &sigmaJ) < 0)
    {
        printf("BondedPair cannot find %s  LJ parameter\n", atmTypeJ);
    }

    *eps = sqrt(epsI * epsJ);
    *sigma = sigmaI + sigmaJ;

    return 0;
}

int genBondPair(CHARMM_PARMS* charmmParms)
{
    int bpairTotalSize = 0;

    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        bpairTotalSize += resiConn->bondListSize;
        bpairTotalSize += resiConn->angleListSize;
        //bpairTotalSize+=resiConn->torsListSize;
    }

    BPAIR_CONN *bpairHeap = ddcMalloc(bpairTotalSize * sizeof (BPAIR_CONN));
    BPAIR_CONN **bpairPtrs = ddcMalloc(bpairTotalSize * sizeof (BPAIR_CONN*));
    for (int i = 0; i < bpairTotalSize; ++i)
    {
        bpairPtrs[i] = &bpairHeap[i];
    }

    CHARMM_PAR* charmmPar = charmmParms->charmmPar;
    bpairTotalSize = 0;

    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];

        resiConn->bpairList = &bpairPtrs[bpairTotalSize];

        for (int i = 0; i < resiConn->bondListSize; ++i)
        {
            bpairPtrs[bpairTotalSize]->atmI = resiConn->bondList[i]->atmI;
            bpairPtrs[bpairTotalSize]->atmJ = resiConn->bondList[i]->atmJ;
            char* atmTypeI = resiConn->bondList[i]->bondPtr->atmTypeI;
            char* atmTypeJ = resiConn->bondList[i]->bondPtr->atmTypeJ;
            double eps;
            double sigma;
            getLJpairParm(charmmPar, atmTypeI, atmTypeJ, &eps, &sigma);
            bpairPtrs[bpairTotalSize]->eps = eps;
            bpairPtrs[bpairTotalSize]->sigma = sigma;

            ++bpairTotalSize;
        }

        // For TIP3 water residue with Shake don't count 1-3 term from angle.
        if (resNameDiffer(resiConn->resName, "TIP3") == 0)
        {
            resiConn->bpairListSize = resiConn->bondListSize;
            continue;
        }

        for (int i = 0; i < resiConn->angleListSize; ++i)
        {
            bpairPtrs[bpairTotalSize]->atmI = resiConn->angleList[i]->atmI;
            bpairPtrs[bpairTotalSize]->atmJ = resiConn->angleList[i]->atmK;
            char* atmTypeI = resiConn->angleList[i]->anglePtr->atmTypeI;
            char* atmTypeJ = resiConn->angleList[i]->anglePtr->atmTypeK;
            double eps;
            double sigma;
            getLJpairParm(charmmPar, atmTypeI, atmTypeJ, &eps, &sigma);
            bpairPtrs[bpairTotalSize]->eps = eps;
            bpairPtrs[bpairTotalSize]->sigma = sigma;

            ++bpairTotalSize;
        }

        /*
                for(int i=0; i<resiConn->torsListSize; ++i){
                    bpairPtrs[bpairTotalSize]->atmI=resiConn->torsList[i]->atmI;
                    bpairPtrs[bpairTotalSize]->atmJ=resiConn->torsList[i]->atmL;
                    char* atmTypeI=resiConn->torsList[i]->atmTypeI;
                    char* atmTypeJ=resiConn->torsList[i]->atmTypeL;
                    double eps;
                    double sigma;
                    getLJpairParm(charmmPar, atmTypeI, atmTypeJ, &eps, &sigma);
                    bpairPtrs[bpairTotalSize]->eps=eps;
                    bpairPtrs[bpairTotalSize]->sigma=sigma;
            
                    ++bpairTotalSize;
                } 
         */

        //resiConn->bpairListSize=resiConn->bondListSize+resiConn->angleListSize+resiConn->torsListSize;
        resiConn->bpairListSize = resiConn->bondListSize + resiConn->angleListSize;

    }

    return 0;
}

int compareBond(const void *v1, const void *v2)
{
    const BOND_CONN *u1 = (BOND_CONN *) v1;
    const BOND_CONN *u2 = (BOND_CONN *) v2;
    if (u1->atmI < u2->atmI)
        return -1;
    else if (u1->atmI > u2->atmI)
        return 1;
    else
        return 0;
}

int compareAngle(const void *v1, const void *v2)
{
    const ANGLE_CONN *u1 = (ANGLE_CONN *) v1;
    const ANGLE_CONN *u2 = (ANGLE_CONN *) v2;
    if (u1->atmJ < u2->atmJ)
        return -1;
    else if (u1->atmJ > u2->atmJ)
        return 1;
    else
        return 0;
}

int compareTors(const void *v1, const void *v2)
{
    const TORS_CONN *u1 = (TORS_CONN *) v1;
    const TORS_CONN *u2 = (TORS_CONN *) v2;
    if (u1->atmJ < u2->atmJ)
        return -1;
    else if (u1->atmJ > u2->atmJ)
        return 1;
    else
        return 0;
}

int compareImpr(const void *v1, const void *v2)
{
    const IMPR_CONN *u1 = (IMPR_CONN *) v1;
    const IMPR_CONN *u2 = (IMPR_CONN *) v2;
    if (u1->atmI < u2->atmI)
        return -1;
    else if (u1->atmI > u2->atmI)
        return 1;
    else
        return 0;
}

int compareBpair(const void *v1, const void *v2)
{
    const BPAIR_CONN *u1 = (BPAIR_CONN *) v1;
    const BPAIR_CONN *u2 = (BPAIR_CONN *) v2;
    if (u1->atmI < u2->atmI)
        return -1;
    else if (u1->atmI > u2->atmI)
        return 1;
    else
        return 0;
}

int compareCmap(const void *v1, const void *v2)
{
    const CMAP_CONN *u1 = (CMAP_CONN *) v1;
    const CMAP_CONN *u2 = (CMAP_CONN *) v2;
    //use N to compare
    if (u1->atmI < u2->atmI)
        return -1;
    else if (u1->atmI > u2->atmI)
        return 1;
    else
        return 0;
}

int compareWeight(const void *v1, const void *v2)
{
    const ZEROWEIGHTS *u1 = (ZEROWEIGHTS *) v1;
    const ZEROWEIGHTS *u2 = (ZEROWEIGHTS *) v2;
    if (u1->atmI < u2->atmI)
        return -1;
    else if (u1->atmI > u2->atmI)
        return 1;
    else
        return 0;
}

void genAtmRange(CHARMM_PARMS* charmmParms)
{
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[i];
        resiConn->atmRanges = ddcMalloc(resiConn->atomListSize * sizeof (ATMRANGE*));
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            ATMRANGE* atmRange = ddcMalloc(sizeof (ATMRANGE));
            resiConn->atmRanges[j] = atmRange;
            // initialize the start indexes
            atmRange->bondRange.start = -1;
            atmRange->angleRange.start = -1;
            atmRange->cosangleRange.start = -1;
            atmRange->rebangleRange.start = -1;
            atmRange->torsRange.start = -1;
            atmRange->imprRange.start = -1;
            atmRange->bpairRange.start = -1;
            atmRange->cmapRange.start = -1;
        }

        BOND_CONN* bondArray = 0;
        ANGLE_CONN* angleArray = 0;
        ANGLE_CONN* cosangleArray = 0;
        ANGLE_CONN* rebangleArray = 0;
        TORS_CONN* torsArray = 0;
        IMPR_CONN* imprArray = 0;
        BPAIR_CONN* bpairArray = 0;
        CMAP_CONN* cmapArray = 0;

        //Solute ions have list sizes = 0
        if (resiConn->bondListSize > 0)
            bondArray = *(resiConn->bondList);
        if (resiConn->angleListSize > 0)
            angleArray = *(resiConn->angleList);
        if (resiConn->cosangleListSize > 0)
            cosangleArray = *(resiConn->cosangleList);        
        if (resiConn->rebangleListSize > 0)
            rebangleArray = *(resiConn->rebangleList);
        if (resiConn->torsListSize > 0)
            torsArray = *(resiConn->torsList);
        if (resiConn->imprListSize > 0)
            imprArray = *(resiConn->imprList);
        if (resiConn->bpairListSize > 0)
            bpairArray = *(resiConn->bpairList);
        if (resiConn->cmapListSize > 0)
            cmapArray = *(resiConn->cmapList);

        if (resiConn->bondListSize > 0)
            qsort(bondArray, resiConn->bondListSize, sizeof (BOND_CONN), compareBond);

        if (resiConn->angleListSize > 0)
            qsort(angleArray, resiConn->angleListSize, sizeof (ANGLE_CONN), compareAngle);

        if (resiConn->cosangleListSize > 0)
            qsort(cosangleArray, resiConn->cosangleListSize, sizeof (ANGLE_CONN), compareAngle);

        if (resiConn->rebangleListSize > 0)
            qsort(rebangleArray, resiConn->rebangleListSize, sizeof (ANGLE_CONN), compareAngle);
        
        if (resiConn->torsListSize > 0)
            qsort(torsArray, resiConn->torsListSize, sizeof (TORS_CONN), compareTors);

        if (resiConn->imprListSize > 0)
            qsort(imprArray, resiConn->imprListSize, sizeof (IMPR_CONN), compareImpr);

        if (resiConn->bpairListSize > 0)
            qsort(bpairArray, resiConn->bpairListSize, sizeof (BPAIR_CONN), compareBpair);

        if (resiConn->cmapListSize > 0)
            qsort(cmapArray, resiConn->cmapListSize, sizeof (CMAP_CONN), compareCmap);

        // Print out sorted array
        /*
        printf("Residue Name: %s  ID: %d  Nter: %d  Cter: %d\n", resiConn->resName, resiConn->resID, resiConn->nTer, resiConn->cTer);
        printf("bond\n");
        for(int j=0; j<resiConn->bondListSize; j++){
            printf("(%d, %d) ", resiConn->bondList[j]->atmI, resiConn->bondList[j]->atmJ);
        }
        printf("\n");
        printf("Angle\n");
        for(int j=0; j<resiConn->angleListSize; j++){
            printf("(%d, %d, %d) ", resiConn->angleList[j]->atmI, resiConn->angleList[j]->atmJ, resiConn->angleList[j]->atmK);
        }
        printf("\n");        
        printf("TORSION\n");
        for(int j=0; j<resiConn->torsListSize; j++){
            printf("(%d, %d, %d, %d) ", resiConn->torsList[j]->atmI, resiConn->torsList[j]->atmJ, resiConn->torsList[j]->atmK, resiConn->torsList[j]->atmL);
        }
        printf("\n");        
         */
        // Assign start and end point for each atom.
        int atomInd = -1;
        for (int j = 0; j < resiConn->bondListSize; j++)
        {
            int atomI = resiConn->bondList[j]->atmI;
            ATMRANGE* atmRange = resiConn->atmRanges[atomI];
            if (atomInd != atomI)
            { // if atomI changes, the range starts here.
                atmRange->bondRange.start = j;
            }
            atmRange->bondRange.end = j + 1; // +1 range ends the place
            atomInd = atomI;
        }

        atomInd = -1;
        for (int j = 0; j < resiConn->angleListSize; j++)
        {
            int atomJ = resiConn->angleList[j]->atmJ;
            ATMRANGE* atmRange = resiConn->atmRanges[atomJ];
            if (atomInd != atomJ)
            { // if atomI changes, the range starts here.
                atmRange->angleRange.start = j;
            }
            atmRange->angleRange.end = j + 1; // +1 range ends the place
            atomInd = atomJ;
        }

        atomInd = -1;
        for (int j = 0; j < resiConn->cosangleListSize; j++)
        {
            int atomJ = resiConn->cosangleList[j]->atmJ;
            ATMRANGE* atmRange = resiConn->atmRanges[atomJ];
            if (atomInd != atomJ)
            { // if atomI changes, the range starts here.
                atmRange->cosangleRange.start = j;
            }
            atmRange->cosangleRange.end = j + 1; // +1 range ends the place
            atomInd = atomJ;
        }
        
        atomInd = -1;
        for (int j = 0; j < resiConn->rebangleListSize; j++)
        {
            int atomJ = resiConn->rebangleList[j]->atmJ;
            ATMRANGE* atmRange = resiConn->atmRanges[atomJ];
            if (atomInd != atomJ)
            { // if atomI changes, the range starts here.
                atmRange->rebangleRange.start = j;
            }
            atmRange->rebangleRange.end = j + 1; // +1 range ends the place
            atomInd = atomJ;
        }
        
        atomInd = -1;
        for (int j = 0; j < resiConn->torsListSize; j++)
        {
            int atomJ = resiConn->torsList[j]->atmJ;
            ATMRANGE* atmRange = resiConn->atmRanges[atomJ];
            if (atomInd != atomJ)
            { // if atomI changes, the range starts here.
                atmRange->torsRange.start = j;
            }
            atmRange->torsRange.end = j + 1; // +1 range ends the place
            atomInd = atomJ;
        }

        atomInd = -1;
        for (int j = 0; j < resiConn->imprListSize; j++)
        {
            int atomI = resiConn->imprList[j]->atmI;
            ATMRANGE* atmRange = resiConn->atmRanges[atomI];
            if (atomInd != atomI)
            { // if atomI changes, the range starts here.
                atmRange->imprRange.start = j;
            }
            atmRange->imprRange.end = j + 1; // +1 range ends the place
            atomInd = atomI;
        }

        atomInd = -1;
        for (int j = 0; j < resiConn->bpairListSize; j++)
        {
            int atomI = resiConn->bpairList[j]->atmI;

            ATMRANGE* atmRange = resiConn->atmRanges[atomI];
            if (atomInd != atomI)
            { // if atomI changes, the range starts here.
                atmRange->bpairRange.start = j;
            }
            atmRange->bpairRange.end = j + 1; // +1 range ends the place
            atomInd = atomI;
        }

        atomInd = -1;
        for (int j = 0; j < resiConn->cmapListSize; j++)
        {
            int atomI = resiConn->cmapList[j]->atmI;
            ATMRANGE* atmRange = resiConn->atmRanges[atomI];
            if (atomInd != atomI)
            {
                atmRange->cmapRange.start = j;
            }
            atmRange->cmapRange.end = j + 1;
            atomInd = atomI;
        }
    }
}

void sortBondList(CHARMM_PARMS*charmmParms)
{
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        if (resiConn->bondListSize == 0)
            continue;
        BOND_CONN* bondArray = *(resiConn->bondList);
        //CHECK: These for loops are <=, changed ddcMalloc
        //int *bondCnt=resiConn->bondCnt = ddcMalloc((resiConn->atomListSize + 1 )*sizeof(int));  // Added 1 element to size to allow for bonds between residues
        int *bondCnt = resiConn->bondCnt = ddcMalloc((resiConn->atomListSize + 2) * sizeof (int)); // Added 1 element to size to allow for bonds between residues
        BNDCNT_LIST *bndCntList = resiConn->bndCntList = ddcMalloc((resiConn->atomListSize + 2) * sizeof (BNDCNT_LIST)); // Tracking bond list for atom.
        //      printf("CHECK: atomListSize = %d\n",resiConn->atomListSize+1);
        for (int ii = 0; ii <= resiConn->atomListSize + 1; ii++)
        {
            bondCnt[ii] = 0;
            for (int jj = 0; jj < ATOMMAXBONDNUM; jj++)
            {
                bndCntList[ii].atom[jj] = -1; //initialize atom bond list
            }
        }

        //CHECK: These for loops are <=, changed ddcMalloc
        //int *firstBond=resiConn->firstBond = ddcMalloc((resiConn->atomListSize + 1 )*sizeof(int)); // Added 1 element so firstBond[atomListSize]=nBonds  can be used to terminal search over bonds 
        int *firstBond = resiConn->firstBond = ddcMalloc((resiConn->atomListSize + 2) * sizeof (int)); // Added 1 element so firstBond[atomListSize]=nBonds  can be used to terminal search over bonds 
        for (int ii = 0; ii <= resiConn->atomListSize + 1; ii++) resiConn->firstBond[ii] = -1;

        for (int ii = 0; ii < resiConn->bondListSize; ii++)
        {
            int atmI = bondArray[ii].atmI;
            int atmJ = bondArray[ii].atmJ;
            bndCntList[atmI].atom[bondCnt[atmI]] = atmJ;
            bndCntList[atmJ].atom[bondCnt[atmJ]] = atmI;
            bondCnt[atmI]++;
            bondCnt[atmJ]++;
            if (atmI > atmJ)
            {
                bondArray[ii].atmI = atmJ;
                bondArray[ii].atmJ = atmI;
            }
        }
        if (resiConn->resType == 1)
        {
            // +N in the next residue has 3 bonds.
            bondCnt[resiConn->atomListSize] = 3;
            bndCntList[resiConn->atomListSize].atom[0] = resiConn->atomListSize - 2; // C
            bndCntList[resiConn->atomListSize].atom[1] = resiConn->atomListSize + 1; //+HN
            bndCntList[resiConn->atomListSize].atom[2] = resiConn->atomListSize + 2; //+CA
            if (resiConn->nTer == 0)
            {
                // if specie name is x or c, current residue N has 3 bonds.
                bondCnt[0] = 3;
                bndCntList[0].atom[0] = -2; // -C
                bndCntList[0].atom[1] = 1; // HN or PRO-CD
                if (resNameDiffer(resiConn->resName, "PRO") == 0)
                {
                    bndCntList[0].atom[2] = 4; // PRO-CA
                }
                else
                {
                    bndCntList[0].atom[2] = 2; // CA
                }
            }
        }

        qsort(bondArray, resiConn->bondListSize, sizeof (BOND_CONN), compareBond);
        int last = -1;
        for (int ii = 0; ii < resiConn->bondListSize; ii++)
        {
            int atmI = bondArray[ii].atmI;
            if (atmI != last)
            {
                firstBond[atmI] = ii;
            }
            last = atmI;
        }
        firstBond[resiConn->atomListSize] = resiConn->bondListSize;
        for (int ii = resiConn->atomListSize - 1; ii > 0; ii--)
        {
            if (firstBond[ii] == -1) firstBond[ii] = firstBond[ii + 1];
        }
    }
}

void findHyrogenHeavy(BOND_CONN** hydrogenbondList, int numHydrogenbond, int atomID, HEAVYATM* heavyAtmPtr)
{

    heavyAtmPtr->numHydrogen = 0;
    for (int i = 0; i < numHydrogenbond; ++i)
    {
        BOND_CONN* bndPtr = hydrogenbondList[i];
        if (bndPtr->atmI == atomID)
        {
            heavyAtmPtr->hydrogen[heavyAtmPtr->numHydrogen] = bndPtr->atmJ;
            heavyAtmPtr->hBndConn[heavyAtmPtr->numHydrogen] = bndPtr;
            heavyAtmPtr->numHydrogen++;
        }
        if (bndPtr->atmJ == atomID)
        {
            heavyAtmPtr->hydrogen[heavyAtmPtr->numHydrogen] = bndPtr->atmI;
            heavyAtmPtr->hBndConn[heavyAtmPtr->numHydrogen] = bndPtr;
            heavyAtmPtr->numHydrogen++;
        }
    }
    if (heavyAtmPtr->numHydrogen > 4)
    {
        printf("foundHyrogenHeavy: find more than 4 hydrogens bonded to %d\n", atomID);
    }
}

int findGrpAtom(RESI_CONN* resiConn, int hydrogenAtmID, int* group, int* grpAtm)
{

    int atomID = 0;
    for (int g = 0; g < resiConn->groupListSize; ++g)
    {
        TGROUP_PARMS *grp = resiConn->groupList[g];

        for (int a = 0; a < grp->grpAtomSize; ++a)
        {
            TATOM_PARMS *atm = grp->grpAtoms[a];

            if (atomID == hydrogenAtmID)
            {
                if (atm->atmTypePtr->mass < HYDROGENMASSCUTOFF)
                {
                    *group = g;
                    *grpAtm = a;
                    return 1;
                }
                else
                {
                    printf("foundGrpAtom: atomID %d is not hydrogen.\n", atomID);
                    return -1;
                }

            }

            atomID++;
        }
    }
    printf("foundGrpAtom: Cannot find atomID %d.\n", atomID);
    return -1;
}

void genHeavyAtms(CHARMM_PARMS* charmmParms)
{

    int totHeavyAtms = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];

        for (int g = 0; g < resiConn->groupListSize; ++g)
        {
            TGROUP_PARMS *grp = resiConn->groupList[g];

            for (int a = 0; a < grp->grpAtomSize; ++a)
            {
                TATOM_PARMS *atm = grp->grpAtoms[a];
                if (atm->atmTypePtr->mass > HYDROGENMASSCUTOFF)
                {
                    totHeavyAtms++;
                }
            }
        }
    }
    // For TIP3 water treat one of the hydrogens as heavy atom
    totHeavyAtms++;

    CHARMM_TOP* charmmTop = charmmParms->charmmTop;

    HEAVYATM *hAtmHeap = ddcMalloc(totHeavyAtms * sizeof (HEAVYATM));
    HEAVYATM **hAtmPtrs = ddcMalloc(totHeavyAtms * sizeof (HEAVYATM*));
    for (int i = 0; i < totHeavyAtms; ++i)
    {
        hAtmPtrs[i] = &hAtmHeap[i];
    }

    int totHeavyAtmCount = 0;
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        int heavyAtmCount = 0;
        resiConn->heavyAtms = &hAtmPtrs[totHeavyAtmCount];

        // If it is a TIP3 water, treat differently.
        // TIP3 water is hard-core assuming the format:
        // RESI TIP3         0.000 ! tip3p water model, generate using noangle nodihedral
        // GROUP
        // ATOM OH2  OT     -0.834
        // ATOM H1   HT      0.417
        // ATOM H2   HT      0.417
        // BOND OH2 H1 OH2 H2 H1 H2    ! the last bond is needed for shake
        if (resNameDiffer(resiConn->resName, "TIP3") == 0)
        {
            //OH2
            HEAVYATM * heavyAtmPtr = hAtmPtrs[totHeavyAtmCount];
            heavyAtmPtr->atomID = 0;
            heavyAtmPtr->atom = 0;
            heavyAtmPtr->group = 0;
            heavyAtmPtr->numHydrogen = 2;
            heavyAtmPtr->hydrogen[0] = 1;
            heavyAtmPtr->hGroup[0] = 0;
            heavyAtmPtr->hGrpAtm[0] = 1;
            heavyAtmPtr->hBndConn[0] = resiConn->bondList[0]; // BOND OH2 H1 
            heavyAtmPtr->hydrogen[1] = 2;
            heavyAtmPtr->hGroup[1] = 0;
            heavyAtmPtr->hGrpAtm[1] = 2;
            heavyAtmPtr->hBndConn[1] = resiConn->bondList[1]; // BOND OH2 H1 
            totHeavyAtmCount++;
            //H1
            heavyAtmPtr = hAtmPtrs[totHeavyAtmCount];
            heavyAtmPtr->atomID = 1;
            heavyAtmPtr->atom = 1;
            heavyAtmPtr->group = 0;
            heavyAtmPtr->numHydrogen = 1;
            heavyAtmPtr->hydrogen[0] = 2;
            heavyAtmPtr->hGroup[0] = 0;
            heavyAtmPtr->hGrpAtm[0] = 2;
            heavyAtmPtr->hBndConn[0] = resiConn->bondList[2]; // BOND H1 H2 
            totHeavyAtmCount++;

            resiConn->heavyAtmsSize = 2;

            continue;

        }

        //Allocate a temporary array for bond containing hydrogen.
        BOND_CONN** hydrogenbondList = ddcMalloc(resiConn->bondListSize * sizeof (BOND_CONN*));
        int numHydrogenbond = 0;
        for (int b = 0; b < resiConn->bondListSize; ++b)
        {
            double massI = getMassbyAtmType(charmmTop, resiConn->bondList[b]->bondPtr->atmTypeI);
            double massJ = getMassbyAtmType(charmmTop, resiConn->bondList[b]->bondPtr->atmTypeJ);
            if (massI < 0 || massJ < 0)
            {
                printf("genHeavyAtms: cannot find mass for bond %s - %s \n",
                       resiConn->bondList[b]->bondPtr->atmTypeI,
                       resiConn->bondList[b]->bondPtr->atmTypeJ);
                continue;
            }
            // 
            if (massI < HYDROGENMASSCUTOFF || massJ < HYDROGENMASSCUTOFF)
            {
                hydrogenbondList[numHydrogenbond] = resiConn->bondList[b];
                numHydrogenbond++;
            }
        }

        int atomID = 0;
        for (int g = 0; g < resiConn->groupListSize; ++g)
        {
            TGROUP_PARMS *grp = resiConn->groupList[g];

            for (int a = 0; a < grp->grpAtomSize; ++a)
            {
                TATOM_PARMS *atm = grp->grpAtoms[a];
                if (atm->atmTypePtr->mass > HYDROGENMASSCUTOFF)
                {
                    HEAVYATM * heavyAtmPtr = hAtmPtrs[totHeavyAtmCount];
                    heavyAtmPtr->atomID = atomID;
                    heavyAtmPtr->group = g;
                    heavyAtmPtr->atom = a;
                    findHyrogenHeavy(hydrogenbondList, numHydrogenbond, atomID, heavyAtmPtr);
                    heavyAtmCount++;
                    totHeavyAtmCount++;
                }

                atomID++;
            }
        }
        resiConn->heavyAtmsSize = heavyAtmCount;

        for (int i = 0; i < resiConn->heavyAtmsSize; i++)
        {
            HEAVYATM * heavyAtmPtr = resiConn->heavyAtms[i];
            for (int h = 0; h < heavyAtmPtr->numHydrogen; h++)
            {
                int group = 0;
                int grpAtm = 0;
                if (findGrpAtom(resiConn, heavyAtmPtr->hydrogen[h], &group, &grpAtm) > 0)
                {
                    heavyAtmPtr->hGroup[h] = group;
                    heavyAtmPtr->hGrpAtm[h] = grpAtm;
                }
                else
                {
                    printf("genHeavyAtms: cannot find group and grpAtm for %s heavy atom %d hydrogen %d.\n",
                           resiConn->resName, heavyAtmPtr->atom, heavyAtmPtr->hydrogen[h]);
                }

            }
        }

        // deallocate the temporary array
        ddcFree(hydrogenbondList);
    }

}

int genConn(CHARMM_PARMS* charmmParms)
{
    genBond(charmmParms);
    sortBondList(charmmParms); // sort bond for quick bond look up
    genAngle(charmmParms);
    genTorsion(charmmParms);
    genBondPair(charmmParms);
    genCmap(charmmParms);
    genWeights(charmmParms);
    genAtmRange(charmmParms); // determine the range of atom in the local task
    genHeavyAtms(charmmParms);

    //numOfBonds=0;
    //numOfAngles=0;
    //numOfDihedral=0;
    //numOfimporper=0;

    return 0;
}

int parseCharmmParms(const char *topFile, const char *parFile, CHARMM_PARMS* charmmParms)
{
    timestamp("Start Bio CHARMM Parameter Parsing");
    charmmParms->charmmPar = (CHARMM_PAR*) ddcMalloc(sizeof (CHARMM_PAR));
    parseCharmmPar(parFile, charmmParms->charmmPar);

    charmmParms->charmmTop = (CHARMM_TOP*) ddcMalloc(sizeof (CHARMM_TOP));
    parseCharmmTop(topFile, charmmParms->charmmTop);

    timestamp("END   Bio CHARMM Parameter Parsing");

    timestamp("Start Bio CHARMM Connectivity Building");

    precomputeCMAP(charmmParms);
    genConn(charmmParms);
    timestamp("END   Bio CHARMM Connectivity Building");
    return 0;
}

RESI_CONN* findResiConnNew(CHARMM_PARMS* charmmParms, char* name)
{
    char resName[16];
    strcpy(resName, name);
    char terminusType[] = {'x', 'n', 'c'};
    char *terminusPtr = strpbrk(resName, "cnx");
    char terminus = *terminusPtr;
    *terminusPtr = '\0';
    for (int r = 0; r < charmmParms->resiConnSize; ++r)
    {
        RESI_CONN* resiConn = charmmParms->resiConnList[r];
        int iTer = resiConn->nTer + 2 * resiConn->cTer;
        if (resNameDiffer(resiConn->resName, resName) == 0 && terminus == terminusType[iTer])
        {
            return resiConn;
        }
    }
    printf("Cannot find residue connectivity for name=%s\n", resName);
    exit(EXIT_FAILURE);
}
