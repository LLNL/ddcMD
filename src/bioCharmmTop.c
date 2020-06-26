#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "ddcMalloc.h"
#include "bioCharmmTop.h"
#include "bioString.h"
#include "units.h"

int resNameDiffer(char* resName, char *name);

typedef enum Top_type
{
    MASS, RESI, GROUP, TOPATOM, TOPBOND, IMPR, CMAP, DELETE, IC, DONOR, ACCEPTOR, PRES, SKIP, TOPEND, SPECIES1, SPECIES2
} TopType;

typedef struct top_size_str
{
    int massSize;
    int resiSize;
    int presSize;
    int grpSize;
    int atmSize;
    int bondSize;
    int imprSize;
    int cmapSize;
    int donorSize;
    int acceptorSize;
    int delatmSize;
    int icSize;
    int species1Size;
    int species2Size;
} TOP_SIZE;

TopType topType(const char* line);
int processMassLine(char* line, MASS_PARMS *massParm);
int processResiLine(char* line, RESI_PARMS *resiParm);
int processAtomLine(char* line, TATOM_PARMS *atomParm, MASS_PARMS **atmMassStack, int massParmSize);
int processBondLine(char* line, TWO_ATOM **bondList, int bNumber);
int processImpOrCmapLine(char* line, FOUR_ATOM **fatomList, int bNumber);
int processDonOrAccLine(char* line, ONE_ATOM **atomList, int bNumber);
int processDelAtmLine(char* line, ONE_ATOM **atomList, int bNumber);
int processICLine(char* line, TIC_PARMS *icParm);
int processSpeciesLine(char* line, ONE_ATOM **atomList, int bNumber);
int setWeightType(CHARMM_TOP* charmTop);

TATOM_PARMS* atmPtrbyName(char* name, TATOM_PARMS **atomList, int atomListSize);

/*
int bondPtr(TWO_ATOM *bondName, TATOM_PARMS **atomList, int atomListSize, TTWO_PARMS *bondParm);
 */

TopType topType(const char* line)
{
    if (line[0] == '\0' || line[0] == ' ' || line[0] == '!') return SKIP;

    char word[9];
    strncpy(word, line, 9);
    word[8] = '\0';
    //MASS, RESI, GROUP, ATOM, IMPR, CMAP, IC, DONOR, ACCEPTOR, PRES
    if (strncmp(word, "MASS", 4) == 0) return MASS;
    else if (strncmp(word, "RESI", 4) == 0) return RESI;
    else if (strncmp(word, "PRES", 4) == 0) return PRES;
    else if (strncmp(word, "GROUP", 5) == 0) return GROUP;
    else if (strncmp(word, "ATOM", 4) == 0) return TOPATOM;
    else if (strncmp(word, "BOND", 4) == 0 || strncmp(word, "DOUBLE", 6) == 0) return TOPBOND;
    else if (strncmp(word, "IMPR", 4) == 0) return IMPR;
    else if (strncmp(word, "CMAP", 4) == 0) return CMAP;
    else if (strncmp(word, "DELETE", 4) == 0) return DELETE;
    else if (strncmp(word, "IC", 2) == 0) return IC;
    else if (strncmp(word, "DONOR", 5) == 0) return DONOR;
    else if (strncmp(word, "ACCEPTOR", 8) == 0) return ACCEPTOR;
    else if (strncmp(word, "SPECIES1", 8) == 0) return SPECIES1;
    else if (strncmp(word, "SPECIES2", 8) == 0) return SPECIES2;
    else if (strncmp(word, "END", 3) == 0) return TOPEND;

    return SKIP;

}

int countWords(const char sentence[ ])
{
    int counted = 0; // result

    // state:
    const char* it = sentence;
    int inword = 0;

    do switch (*it)
        {
        case '\0':
        case ' ': case '\t': case '\n': case '\r': // TODO others?
            if (inword)
            {
                inword = 0;
                counted++;
            }
            break;
        default: inword = 1;
        }
    while (*it++);

    return counted;
}
static int massParmNum = 0;
static int resiParmNum = 0;

static int grpNum = 0;
static int totalGrpNum = 0;
static int atmNum = 0;
static int totalAtmNum = 0;
static int grpAtmNum = 0;
static int bondNum = 0;
static int totalBondNum = 0;
static int imprNum = 0;
static int totalImprNum = 0;
static int cmapNum = 0;
//Note: total Cmap Num refers to the total number of dihedrals
//the actual number of Cmap calculations is totalCmapNum/2
static int totalCmapNum = 0;
static int donNum = 0;
static int totalDonNum = 0;
static int accNum = 0;
static int totalAccNum = 0;
static int delNum = 0;
static int totalDelNum = 0;
static int icNum = 0;
static int totalICNum = 0;

static int spec1Num = 0;
static int totalSpec1Num = 0;
static int spec2Num = 0;
static int totalSpec2Num = 0;

static int firstRes = 0;
static int firstGrp = 0;

static int firstResAtm = 0;
static int firstGrpAtm = 0;

static int firstBond = 0;
static int firstImpr = 0;
static int firstCmap = 0;
static int firstDon = 0;
static int firstAcc = 0;
static int firstDel = 0;
static int firstIC = 0;

static int firstSpec1 = 0;
static int firstSpec2 = 0;

void resiReset(CHARMM_TOP *charmmTops, TGROUP_PARMS **groupPtrs)
{
    RESI_PARMS *prevResi = charmmTops->resiParms[resiParmNum - 1];
    // Save the previous residue infos from stacks.
    //GROUP
    if (grpNum > 0)
    {
        TGROUP_PARMS *prevGrp = groupPtrs[totalGrpNum - 1];
        prevGrp->grpAtomSize = grpAtmNum;
        grpAtmNum = 0;
    }
    prevResi->groupListSize = grpNum;
    grpNum = 0;
    firstGrp = 0;
    prevResi->atomListSize = atmNum; //ATOM
    atmNum = 0;
    prevResi->bondListSize = bondNum; //BOND
    bondNum = 0;
    firstBond = 0;
    prevResi->imprListSize = imprNum; //IMPR
    imprNum = 0;
    firstImpr = 0;
    prevResi->cmapListSize = cmapNum; //CMAP
    cmapNum = 0;
    firstCmap = 0;
    prevResi->donorListSize = donNum; //DONOR
    donNum = 0;
    firstDon = 0;
    prevResi->accepListSize = accNum; //ACCEPTOR
    accNum = 0;
    firstAcc = 0;
    prevResi->delAtmListSize = delNum; //DELETE ATOM
    delNum = 0;
    firstDel = 0;
    prevResi->icListSize = icNum; //IC
    icNum = 0;
    firstIC = 0;
    prevResi->species1Size = spec1Num; //Transforming Species 1
    spec1Num = 0;
    firstSpec1 = 0;
    prevResi->species2Size = spec2Num; //Transforming Species 2
    spec2Num = 0;
    firstSpec2 = 0;
    //END
}

int getTopSize(const char *fileName, TOP_SIZE *topSize)
{

    int massParmSize = 0;
    int resiParmSize = 0;
    int presParmSize = 0;

    int grpSize = 0;
    int atmSize = 0;
    int bondSize = 0;
    int imprSize = 0;
    //Number of cmap calculations = cmapSize/2  
    int cmapSize = 0;
    int donorSize = 0;
    int acceptorSize = 0;
    int delatmSize = 0;
    int icSize = 0;

    int species1Size = 0;
    int species2Size = 0;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1)
    {
        //printf("Retrieved line of length %zu :\n", read);
        TopType lt = topType(line);
        if (lt == SKIP)
        {
            continue;
        }

        char *noCLine = stripComment(line);
        line = rstrip(noCLine);
        //printf("Line: %s \n", line);

        if (lt == MASS)
        {
            massParmSize++;
        }
        else if (lt == RESI)
        {
            resiParmSize++;
        }
        else if (lt == PRES)
        {
            presParmSize++;
        }
        else if (lt == GROUP)
        {
            grpSize++;
        }
        else if (lt == TOPATOM)
        {
            atmSize++;
        }
        else if (lt == TOPBOND)
        {
            int bondResSize = countWords(line);
            if (bondResSize % 2 != 1) printf("Bond has odd number of atoms: %s \n", line);
            bondSize = bondSize + (bondResSize - 1) / 2;
        }
        else if (lt == IMPR)
        {
            int imprResSize = countWords(line);
            if (imprResSize % 4 != 1) printf("IMPR does not have 4 atoms: %s \n", line);
            imprSize = imprSize + (imprResSize - 1) / 4;
        }
        else if (lt == CMAP)
        {
            int cmapResSize = countWords(line);
            cmapSize = cmapSize + (cmapResSize - 1) / 4;
        }
        else if (lt == DONOR)
        {
            donorSize = donorSize + countWords(line) - 1;
        }
        else if (lt == ACCEPTOR)
        {
            acceptorSize = acceptorSize + countWords(line) - 1;
        }
        else if (lt == IC)
        {
            icSize++;
        }
        else if (lt == DELETE)
        {
            delatmSize++;
        }
        else if (lt == SPECIES1)
        {
            species1Size++;
        }
        else if (lt == SPECIES2)
        {
            species2Size++;
        }
    }

    topSize->massSize = massParmSize;
    topSize->resiSize = resiParmSize;
    topSize->presSize = presParmSize;
    topSize->grpSize = grpSize;
    topSize->atmSize = atmSize;
    topSize->bondSize = bondSize;
    topSize->imprSize = imprSize;
    topSize->cmapSize = cmapSize;
    topSize->donorSize = donorSize;
    topSize->acceptorSize = acceptorSize;
    topSize->icSize = icSize;
    topSize->delatmSize = delatmSize;

    topSize->species1Size = species1Size;
    topSize->species2Size = species2Size;

    printf("CHECK: counted Species1 = %d Species2 = %d\n", species1Size, species2Size);

    fclose(fp);
    if (line)
        free(line);

    return 0;
}

int parseCharmmTop(const char *fileName, CHARMM_TOP *charmmTops)
{

    TOP_SIZE topSize;
    getTopSize(fileName, &topSize);

    MASS_PARMS *atmMassStack = ddcMalloc(topSize.massSize * sizeof (MASS_PARMS));
    charmmTops->massParmSize = topSize.massSize;
    charmmTops->massParms = ddcMalloc(topSize.massSize * sizeof (MASS_PARMS*));
    for (int i = 0; i < topSize.massSize; ++i)
    {
        charmmTops->massParms[i] = &atmMassStack[i];
    }

    int resipresSize = topSize.resiSize + topSize.presSize;
    RESI_PARMS *resParmStack = ddcMalloc(resipresSize * sizeof (RESI_PARMS));
    RESI_PARMS **resParmPtrs = ddcMalloc(resipresSize * sizeof (RESI_PARMS*));
    //RESI
    charmmTops->resiParmSize = topSize.resiSize;
    charmmTops->resiParms = &resParmPtrs[0];
    for (int i = 0; i < topSize.resiSize; ++i)
    {
        charmmTops->resiParms[i] = &resParmStack[i];
    }
    //PRES
    charmmTops->presParmSize = topSize.presSize;
    charmmTops->presParms = &resParmPtrs[topSize.resiSize];
    for (int i = 0; i < topSize.presSize; ++i)
    {
        charmmTops->presParms[i] = &resParmStack[topSize.resiSize + i];
    }

    TGROUP_PARMS *groupStack = ddcMalloc(topSize.grpSize * sizeof (TGROUP_PARMS));
    TGROUP_PARMS **groupPtrs = ddcMalloc(topSize.grpSize * sizeof (TGROUP_PARMS*));
    for (int i = 0; i < topSize.grpSize; ++i)
    {
        groupPtrs[i] = &groupStack[i];
    }

    //TATOM_PARMS *grpAtomStack=ddcMalloc(topSize.massSize*sizeof(MASS_PARMS));
    TATOM_PARMS *atmParmStack = ddcMalloc(topSize.atmSize * sizeof (TATOM_PARMS));
    TATOM_PARMS **atmParmPtrs = ddcMalloc(topSize.atmSize * sizeof (TATOM_PARMS*));
    for (int i = 0; i < topSize.atmSize; ++i)
    {
        atmParmPtrs[i] = &atmParmStack[i];
    }

    TWO_ATOM *bondStack = ddcMalloc(topSize.bondSize * sizeof (TWO_ATOM));
    TWO_ATOM **bondPtrs = ddcMalloc(topSize.bondSize * sizeof (TWO_ATOM*));
    for (int i = 0; i < topSize.bondSize; ++i)
    {
        bondPtrs[i] = &bondStack[i];
    }

    FOUR_ATOM *imprStack = ddcMalloc(topSize.imprSize * sizeof (FOUR_ATOM));
    FOUR_ATOM **imprPtrs = ddcMalloc(topSize.imprSize * sizeof (FOUR_ATOM*));
    for (int i = 0; i < topSize.imprSize; ++i)
    {
        imprPtrs[i] = &imprStack[i];
    }

    FOUR_ATOM *cmapStack = ddcMalloc(topSize.cmapSize * sizeof (FOUR_ATOM));
    FOUR_ATOM **cmapPtrs = ddcMalloc(topSize.cmapSize * sizeof (FOUR_ATOM*));
    for (int i = 0; i < topSize.cmapSize; ++i)
    {
        cmapPtrs[i] = &cmapStack[i];
    }

    ONE_ATOM *donStack = ddcMalloc(topSize.donorSize * sizeof (ONE_ATOM)); // Convert to atom ptr later
    ONE_ATOM **donPtrs = ddcMalloc(topSize.donorSize * sizeof (ONE_ATOM*));
    for (int i = 0; i < topSize.donorSize; ++i)
    {
        donPtrs[i] = &donStack[i];
    }

    ONE_ATOM *accStack = ddcMalloc(topSize.acceptorSize * sizeof (ONE_ATOM)); // Convert to atom ptr later
    ONE_ATOM **accPtrs = ddcMalloc(topSize.acceptorSize * sizeof (ONE_ATOM*));
    for (int i = 0; i < topSize.acceptorSize; ++i)
    {
        accPtrs[i] = &accStack[i];
    }

    ONE_ATOM *delAtmStack = ddcMalloc(topSize.delatmSize * sizeof (ONE_ATOM)); // Convert to atom ptr later
    ONE_ATOM **delAtmPtrs = ddcMalloc(topSize.delatmSize * sizeof (ONE_ATOM*));
    for (int i = 0; i < topSize.delatmSize; ++i)
    {
        delAtmPtrs[i] = &delAtmStack[i];
    }

    ONE_ATOM *spec1Stack = ddcMalloc(topSize.species1Size * sizeof (ONE_ATOM)); // Convert to atom ptr later
    ONE_ATOM **spec1Ptrs = ddcMalloc(topSize.species1Size * sizeof (ONE_ATOM*));
    for (int i = 0; i < topSize.species1Size; ++i)
    {
        spec1Ptrs[i] = &spec1Stack[i];
    }

    ONE_ATOM *spec2Stack = ddcMalloc(topSize.species2Size * sizeof (ONE_ATOM)); // Convert to atom ptr later
    ONE_ATOM **spec2Ptrs = ddcMalloc(topSize.species2Size * sizeof (ONE_ATOM*));
    for (int i = 0; i < topSize.species2Size; ++i)
    {
        spec2Ptrs[i] = &spec2Stack[i];
    }

    TIC_PARMS *icStack = ddcMalloc(topSize.icSize * sizeof (TIC_PARMS));
    TIC_PARMS **icPtrs = ddcMalloc(topSize.icSize * sizeof (TIC_PARMS*));
    FOUR_ATOM *icAtmsStack = ddcMalloc(topSize.icSize * sizeof (FOUR_ATOM));
    for (int i = 0; i < topSize.icSize; ++i)
    {
        icPtrs[i] = &icStack[i];
        icPtrs[i]->atmsPtr = &icAtmsStack[i];
    }



    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1)
    {
        TopType lt = topType(line);

        if (lt == SKIP) continue;
        char *noCLine = stripComment(line);
        line = rstrip(noCLine);
        switch (lt)
        {
        case MASS:
        {
            MASS_PARMS *massParm = charmmTops->massParms[massParmNum];
            if (processMassLine(line, massParm) < 0) printf("Mass error: %s\n", line);
            massParmNum++;
        }
            break;
        case TOPEND:
            resiReset(charmmTops, groupPtrs);
            break;
        case RESI:
        case PRES:
        {
            if (firstRes != 0) resiReset(charmmTops, groupPtrs);
            RESI_PARMS *resiParm = charmmTops->resiParms[resiParmNum];
            if (processResiLine(line, resiParm) < 0)
            {
                printf("Residue error: %s\n", line);
            }
            if (resNameDiffer(resiParm->resName, "NTER") == 0)
            {
                printf("Debug here resName=: %s\n", resiParm->resName); //resClean
            }
            resiParm->resID = resiParmNum;
            resiParmNum++;
            firstRes = 1;
            firstResAtm = 0;
        }
            break;
        case GROUP:
        {
            if (firstGrp != 0)
            {
                TGROUP_PARMS *prevGrp = groupPtrs[totalGrpNum - 1];
                prevGrp->grpAtomSize = grpAtmNum;
                grpAtmNum = 0;
            }
            TGROUP_PARMS *grpParm = groupPtrs[totalGrpNum];
            grpParm->grpID = grpNum;

            if (firstGrp == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->groupList = &groupPtrs[totalGrpNum];
                firstGrp = 1;
            }

            totalGrpNum++;
            grpNum++;
            firstGrpAtm = 0;
        }
            break;
        case TOPATOM:
        {
            TATOM_PARMS *atomParm = atmParmPtrs[totalAtmNum];
            if (processAtomLine(line, atomParm, charmmTops->massParms, charmmTops->massParmSize) < 0)
            {
                printf("ATOM error: %s", line);
            }

            if (firstGrpAtm == 0)
            {
                groupPtrs[totalGrpNum - 1]->grpAtoms = &atmParmPtrs[totalAtmNum];
                firstGrpAtm = 1;
            }

            if (firstResAtm == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->atomList = &atmParmPtrs[totalAtmNum];
                firstResAtm = 1;
            }

            atmNum++;
            totalAtmNum++;

            atomParm->atmID = grpAtmNum;
            grpAtmNum++;
        }
            break;
        case TOPBOND:
        {
            TWO_ATOM **bondTempPtr = &bondPtrs[totalBondNum];
            if (firstBond == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->bondList = bondTempPtr;
                firstBond = 1;
            }
            int bondTempSize = (countWords(line) - 1) / 2;
            if (processBondLine(line, bondTempPtr, bondTempSize) < 0)
            {
                printf("BOND error: %s\n", line);
            }
            bondNum = bondNum + bondTempSize;
            totalBondNum = totalBondNum + bondTempSize;
        }
            break;
        case IMPR:
        {
            FOUR_ATOM **imprTempPtr = &imprPtrs[totalImprNum];
            if (firstImpr == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->imprList = imprTempPtr;
                firstImpr = 1;
            }
            int imprTempSize = (countWords(line) - 1) / 4;
            if (processImpOrCmapLine(line, imprTempPtr, imprTempSize) < 0)
            {
                printf("Impr error: %s\n", line);
            }
            imprNum = imprNum + imprTempSize;
            totalImprNum = totalImprNum + imprTempSize;
        }
            break;
        case CMAP:
        {
            FOUR_ATOM **cmapTempPtr = &cmapPtrs[totalCmapNum];
            if (firstCmap == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->cmapList = cmapTempPtr;
                firstCmap = 1;
            }
            int cmapTempSize = (countWords(line) - 1) / 4;
            if (processImpOrCmapLine(line, cmapTempPtr, cmapTempSize) < 0)
            {
                printf("CMap error: %s\n", line);
            }
            cmapNum = cmapNum + cmapTempSize;
            totalCmapNum = totalCmapNum + cmapTempSize;
        }
            break;
        case DONOR:
        {
            ONE_ATOM** donTempPtr = &donPtrs[totalDonNum];
            if (firstDon == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->donorList = donTempPtr;
                firstDon = 1;
            }
            int donTempSize = countWords(line) - 1;
            processDonOrAccLine(line, donTempPtr, donTempSize);
            donNum = donNum + donTempSize;
            totalDonNum = totalDonNum + donTempSize;
        }
            break;
        case ACCEPTOR:
        {
            ONE_ATOM** accTempPtr = &accPtrs[totalAccNum];
            if (firstAcc == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->accepList = accTempPtr;
                firstAcc = 1;
            }
            int accTempSize = countWords(line) - 1;
            processDonOrAccLine(line, accTempPtr, accTempSize);
            accNum = accNum + accTempSize;
            totalAccNum = totalAccNum + accTempSize;
        }
            break;
        case DELETE:
        {
            ONE_ATOM** delTempPtr = &delAtmPtrs[totalDelNum];
            if (firstDel == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->delAtmList = delTempPtr;
                firstDel = 1;
            }
            int delTempSize = countWords(line) - 2;
            processDelAtmLine(line, delTempPtr, delTempSize);
            delNum = delNum + delTempSize;
            totalDelNum = totalDelNum + delTempSize;
        }
            break;
        case IC:
        {
            TIC_PARMS *icParm = icPtrs[totalICNum];
            if (firstIC == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->icList = &icPtrs[totalICNum];
                firstIC = 1;
            }
            if (processICLine(line, icParm) < 0)
            {
                printf("IC error: %s", line);
            }
            icNum++;
            totalICNum++;
        }
            break;
        case SPECIES1:
        {
            ONE_ATOM** spec1TempPtr = &spec1Ptrs[totalSpec1Num];
            if (firstSpec1 == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->species1List = spec1TempPtr;
                firstSpec1 = 1;
            }
            int spec1TempSize = countWords(line) - 1;
            printf("CHECK: Found temp Species 1 Size = %d\n", spec1TempSize);
            processSpeciesLine(line, spec1TempPtr, spec1TempSize);
            spec1Num += spec1TempSize;
            totalSpec1Num += spec1TempSize;
            printf("CHECK: total Species 1 Num = %d\n", totalSpec1Num);
        }
            break;
        case SPECIES2:
        {
            ONE_ATOM** spec2TempPtr = &spec2Ptrs[totalSpec2Num];
            if (firstSpec2 == 0)
            {
                charmmTops->resiParms[resiParmNum - 1]->species2List = spec2TempPtr;
                firstSpec2 = 1;
            }
            int spec2TempSize = countWords(line) - 1;
            printf("CHECK: Found temp Species 2 Size = %d\n", spec2TempSize);
            processSpeciesLine(line, spec2TempPtr, spec2TempSize);
            spec2Num += spec2TempSize;
            totalSpec2Num += spec2TempSize;
        }
            break;
        default:
        {
        }
        }
    }

    fclose(fp);
    if (line)
        free(line);

    setWeightType(charmmTops);

    return 0;
}

int processMassLine(char* line, MASS_PARMS *massParm)
{
    //  MASS     1 H      1.00800 H ! polar H    --- CHARMM22
    //  MASS  -1  H          1.00800 ! polar H   --- CHARMM36

    char *token;
    token = strtok(line, " ");
    int count = 0;
    /* walk through other tokens */
    while (token != NULL)
    {
        //printf("\t %d: %s\n", count, token);
        if (count == 2)
        {
            strncpy(massParm->atmType, token, 4);
            massParm->atmType[4] = '\0';
        }
        if (count == 3)
        {
            double mass;
            sscanf(token, "%lf", &mass);
            massParm->mass = units_convert(mass, "m", NULL); // Convert mass unit to default unit.
        }
        if (count == 4)
        {
            strncpy(massParm->element, token, 2);
            massParm->element[2] = '\0';
            return 0;
        }
        token = strtok(NULL, " ");
        count++;
    }

    return 0;
    /*
       int numscan=-1;
       int id=0;

       numscan = sscanf(line, "MASS%7d%4s %10lf%2s",
       &id, massParm->atmType, &massParm->mass, massParm->element);    
     */

}


#if 0 

int processResiLine(char* line, RESI_PARMS *resiParm)
{

    char *stdAA[] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

    //RESI ALA          0.00
    char *token;
    token = strtok(line, " ");
    int count = 0;
    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 1)
        {
            strncpy(resiParm->resName, token, 4);
            resiParm->resName[4] = '\0';
            resiParm->resType = 0;
            int i;
            for (i = 0; i < 22; ++i)
            {
                if (resNameDiffer(resiParm->resName, stdAA[i]) == 0)
                {
                    resiParm->resType = 1;
                    break;
                }
            }
        }
        if (count == 2)
        {
            double charge;
            sscanf(token, "%lf", &charge);
            resiParm->charge = units_convert(charge, "i*t", NULL);
            return 0;
        }
        token = strtok(NULL, " ");
        count++;
    }
    return -1;
}
#endif

int processResiLine(char* line, RESI_PARMS *resiParm)
{

    char *stdAA[] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

    //RESI ALA          0.00
    char *token;
    token = strtok(line, " ");
    if (token == NULL) return -1;
    token = strtok(NULL, " ");
    if (token == NULL) return -1;
    strcpy(resiParm->resName, token);
    resiParm->resType = 0;
    for (int i = 0; i < 22; ++i)
        if (resNameDiffer(resiParm->resName, stdAA[i]) == 0)
        {
            resiParm->resType = 1;
            break;
        }
    token = strtok(NULL, " ");
    if (token == NULL) return -1;
    resiParm->charge *= atof(token) * units_convert(1.0, "i*t", NULL);
    return 0;
}

int processAtomLine(char* line, TATOM_PARMS *atomParm, MASS_PARMS **atmMassStack, int massParmSize)
{
    //ATOM N    NH1    -0.47  !     |
    char *token;
    token = strtok(line, " ");
    int count = 0;

    char atmType[5];
    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 1)
        {
            strncpy(atomParm->atmName, token, 4);
            atomParm->atmName[4] = '\0';
        }
        if (count == 2)
        {
            strncpy(atmType, token, 4);
            atmType[4] = '\0';
            int foundAtmType = 0;
            for (int i = 0; i < massParmSize; i++)
            {
                if (strncmp(atmMassStack[i]->atmType, atmType, 4) == 0)
                {
                    atomParm->atmTypePtr = atmMassStack[i];
                    foundAtmType = 1;
                    break;
                }
            }
            if (foundAtmType == 0)
            {
                printf("Cannot find atom type pinter for: %s", atmType);
            }
        }
        if (count == 3)
        {
            double charge;
            sscanf(token, "%lf", &charge);
            atomParm->charge = units_convert(charge, "i*t", NULL);
            return 0;
        }
        token = strtok(NULL, " ");
        count++;
    }
    return -1;
}

int processBondLine(char* line, TWO_ATOM **bondList, int bNumber)
{
    // BOND CB CA  CG CB   ND2 CG 
    /*
       char *p=line;
       char *newLine=stripComment((p+5));
     */
    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);

    char *token;
    token = strtok(newLine, " ");
    int count = 0;

    ONE_ATOM atomStack[20];
    /* walk through other tokens */
    while (token != NULL)
    {
        strncpy(atomStack[count].atomI, token, 4);
        atomStack[count].atomI[4] = '\0';
        token = strtok(NULL, " ");
        count++;
    }
    if ((count - 1) % 2 != 0)
    {
        printf("BOND pairs are not even, count-1=%d\n", count - 1);
        return -1;
    }
    if (bNumber != ((count - 1) / 2))
    {
        printf("BOND pairs are not equal.\n");
        return -1;
    }
    int i;
    for (i = 0; i < bNumber; i++)
    {
        TWO_ATOM *atomParm = bondList[i];
        strcpy(atomParm->atomI, atomStack[2 * i + 1].atomI);
        strcpy(atomParm->atomJ, atomStack[2 * i + 2].atomI);
    }

    return 0;
}

int processImpOrCmapLine(char* line, FOUR_ATOM **fatomList, int bNumber)
{

    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);
    char *token;
    token = strtok(newLine, " ");
    int count = 0;

    ONE_ATOM atomStack[20];
    /* walk through other tokens */
    while (token != NULL)
    {
        strncpy(atomStack[count].atomI, token, 4);
        atomStack[count].atomI[4] = '\0';
        token = strtok(NULL, " ");
        count++;
    }
    if ((count - 1) % 4 != 0)
    {
        printf("Impr or Cmap list are not four, count-2=%d\n", count - 2);
        return -1;
    }

    if (((count - 1) / 4) != bNumber)
    {
        printf("Impr or Cmap list are not equal.\n");
        return -1;
    }

    for (int i = 0; i < bNumber; i++)
    {
        FOUR_ATOM *atomParm = fatomList[i];
        strcpy(atomParm->atomI, atomStack[4 * i + 1].atomI);
        strcpy(atomParm->atomJ, atomStack[4 * i + 2].atomI);
        strcpy(atomParm->atomK, atomStack[4 * i + 3].atomI);
        strcpy(atomParm->atomL, atomStack[4 * i + 4].atomI);
    }

    return 0;
}

int processDonOrAccLine(char* line, ONE_ATOM **atomList, int bNumber)
{
    //DONOR HE NE
    //DONOR HH11 NH1

    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);

    char *token;
    token = strtok(newLine, " ");
    int count = 0;

    /* walk through other tokens */
    while (token != NULL)
    {
        if (count > 0)
        {
            strncpy(atomList[count - 1]->atomI, token, 4);
            atomList[count - 1]->atomI[4] = '\0';
        }
        token = strtok(NULL, " ");
        count++;
    }
    if (bNumber != (count - 1))
    {
        printf("Donor or Acceptor atom lists are not equal.\n");
    }

    return 0;
}

int processDelAtmLine(char* line, ONE_ATOM **atomList, int bNumber)
{
    //DELETE ATOM HN

    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);

    char *token;
    token = strtok(newLine, " ");
    int count = 0;

    /* walk through other tokens */
    while (token != NULL)
    {
        if (count > 1)
        {
            strncpy(atomList[count - 2]->atomI, token, 4);
            atomList[count - 2]->atomI[4] = '\0';
            break;
        }
        token = strtok(NULL, " ");
        count++;
    }

    if (bNumber != (count - 1))
    {
        printf("Delete atom lists are not equal.\n");
    }

    return 0;
}

int processICLine(char* line, TIC_PARMS *icParm)
{
    //IC -C   CA   *N   HN    1.3496 122.4500  180.0000 116.6700  0.9973        
    FOUR_ATOM *atoms = icParm->atmsPtr;

    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);

    char *token;
    token = strtok(newLine, " ");
    int count = 0;

    // TODO: CHARMM IC Table units are not converted yet.
    while (token != NULL)
    {
        if (count == 1)
        {
            strncpy(atoms->atomI, token, 4);
            atoms->atomI[4] = '\0';
        }
        if (count == 2)
        {
            strncpy(atoms->atomJ, token, 4);
            atoms->atomJ[4] = '\0';
        }
        if (count == 3)
        {
            strncpy(atoms->atomK, token, 4);
            atoms->atomK[4] = '\0';
        }
        if (count == 4)
        {
            strncpy(atoms->atomL, token, 4);
            atoms->atomL[4] = '\0';
        }
        if (count == 5)
        {
            sscanf(token, "%lf", &icParm->kconst1);
        }
        if (count == 6)
        {
            sscanf(token, "%lf", &icParm->angle1);
        }
        if (count == 7)
        {
            sscanf(token, "%lf", &icParm->torsion);
        }
        if (count == 8)
        {
            sscanf(token, "%lf", &icParm->angle2);
        }
        if (count == 9)
        {
            sscanf(token, "%lf", &icParm->kconst2);
            return 0;
        }
        token = strtok(NULL, " ");
        count++;
    }
    return -1;
}

int processSpeciesLine(char* line, ONE_ATOM **atomList, int size)
{

    char *noCLine = stripComment(line);
    char *newLine = rstrip(noCLine);

    char *token;
    token = strtok(newLine, " ");
    int count = 0;
    while (token != NULL)
    {
        if (count > 0)
        {
            strncpy(atomList[count - 1]->atomI, token, 4);
            atomList[count - 1]->atomI[4] = '\0';
            printf("CHECK: species line: %s\n", atomList[count - 1]->atomI);
        }
        token = strtok(NULL, " ");
        count++;
    }

    return 0;
}

int setWeightType(CHARMM_TOP* charmmTop)
{
    for (int r = 0; r < charmmTop->resiParmSize; ++r)
    {
        RESI_PARMS* parms = charmmTop->resiParms[r];
        //printf("CHECK: species1Size = %d species2Size=%d\n", parms->species1Size, parms->species2Size); 

        for (int i = 0; i < parms->atomListSize; i++)
        {

            parms->atomList[i]->species = 0;

            if (parms->species1Size != 0)
            {
                char *atmName = parms->atomList[i]->atmName;
                for (int j = 0; j < parms->species1Size; j++)
                {
                    if (strcmp(atmName, parms->species1List[j]->atomI) == 0)
                    {
                        parms->atomList[i]->species = 1;
                        continue;
                    }
                }
            }
            if (parms->species2Size != 0 && parms->atomList[i]->species == 0)
            {
                char * atmName = parms->atomList[i]->atmName;
                for (int j = 0; j < parms->species2Size; j++)
                {
                    if (strcmp(atmName, parms->species2List[j]->atomI) == 0)
                    {
                        parms->atomList[i]->species = 2;
                        continue;
                    }
                }
            }
            //printf("CHECK: residue = %s species = %d\n", parms->resName, parms->atomList[i]->species);    /resClean
        }
    }
    return 0;
}

/*
   int bondPtr(TWO_ATOM *bondName, TATOM_PARMS **atomList, int atomListSize, TTWO_PARMS *bondParm){
   if(atmPtrbyName(bondName->atomI, atomList, atomListSize, bondParm->atomI)<0){
   printf("BOND AtomI error\n");
   return -1;
   }
   if(atmPtrbyName(bondName->atomJ, atomList, atomListSize, bondParm->atomJ)<0){
   printf("BOND AtomJ error\n");
   return -1;
   }    
   return 0;
   }
 */

TATOM_PARMS* atmPtrbyName(char* name, TATOM_PARMS **atomList, int atomListSize)
{

    for (int i = 0; i < atomListSize; i++)
    {
        if (strcmp(atomList[i]->atmName, name) == 0)
        {
            return atomList[i];
        }
    }
    return 0;
}

int freeCharmmTop(CHARMM_TOP *charmmTops)
{

    for (int i = 0; i < charmmTops->massParmSize; i++)
    {
        free(charmmTops->massParms[i]);
    }
    free(charmmTops->massParms);

    for (int i = 0; i < charmmTops->resiParmSize; i++)
    {
        RESI_PARMS *resiParm = charmmTops->resiParms[i];

        for (int j = 0; j < resiParm->groupListSize; j++)
        {
            TGROUP_PARMS* groupParm = resiParm->groupList[j];

            for (int k = 0; k < groupParm->grpAtomSize; k++)
            {
                free(groupParm->grpAtoms[k]);
            }
            free(groupParm);
        }
        free(resiParm->groupList);

        for (int j = 0; j < resiParm->bondListSize; j++)
        {
            free(resiParm->bondList[j]);
        }
        free(resiParm->bondList);

        for (int j = 0; j < resiParm->imprListSize; j++)
        {
            free(resiParm->imprList[j]);
        }
        free(resiParm->imprList);

        for (int j = 0; j < resiParm->cmapListSize; j++)
        {
            free(resiParm->cmapList[j]);
        }
        free(resiParm->cmapList);

        for (int j = 0; j < resiParm->icListSize; j++)
        {
            free(resiParm->icList[j]->atmsPtr);
            free(resiParm->icList[j]);
        }
        free(resiParm->icList);

        free(resiParm->donorList);
        free(resiParm->accepList);

        free(resiParm->species1List);
        free(resiParm->species2List);

        free(resiParm);
    }
    free(charmmTops->resiParms);

    return 0;
}

RESI_PARMS* findResiParms(CHARMM_TOP* charmmTop, char* resName) //resClean
{

    for (int i = 0; i < charmmTop->resiParmSize; ++i)
    {
        RESI_PARMS *resiParms = charmmTop->resiParms[i];
        if (resNameDiffer(resiParms->resName, resName) == 0)
        {
            return resiParms;
        }
    }

    for (int i = 0; i < charmmTop->presParmSize; ++i)
    {
        RESI_PARMS *resiParms = charmmTop->presParms[i];
        if (resNameDiffer(resiParms->resName, resName) == 0)
        {
            return resiParms;
        }
    }

    return 0;

}

double getMassbyAtmType(CHARMM_TOP* charmmTop, char* atmType)
{
    for (int i = 0; i < charmmTop->massParmSize; ++i)
    {
        MASS_PARMS *massParms = charmmTop->massParms[i];
        if (strcmp(massParms->atmType, atmType) == 0)
        {
            return massParms->mass;
        }
    }
    return -1;
}
