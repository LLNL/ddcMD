#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "ddcMalloc.h"
#include "bioCharmmPar.h"
#include "bioString.h"
#include "codata.h"
#include "units.h"

typedef enum Par_type
{
    BONDPAR, ANGLEPAR, TORSIONPAR, IMPRPAR, CMAPPAR, LJPAR, HBONDPAR, SKIPPAR, CONTENTPAR, ENDPAR
} ParType;

int processBondContent(char* line, BOND_PARMS *bondParm);
int processAngleContent(char* line, ANGLE_PARMS *angleParm);
int processTorContent(char* line, TORSION_PARMS *torParm);
int processImpContent(char* line, IMPROPER_PARMS *imprParm);
int processLJContent(char* line, LJCH_PARMS *ljParm);
int processCMAPContent(char* line, CMAP_PARMS *cmapParm, int gridStart, int* nGridTerms);
int printCMAPContent(CMAP_PARMS *cmapParm);
int getParSize(const char *fileName, int* maxbondptr, int* maxangleptr, int* maxtorptr, int* maximpptr, int* maxljptr, int *maxcmapptr);
int skipLine(const char *line);

ParType parType(const char* line)
{

    char word[9];
    strncpy(word, line, 9);
    word[8] = '\0';
    //MASS, RESI, GROUP, ATOM, IMPR, CMAP, IC, DONOR, ACCEPTOR, PRES
    if (skipLine(line)) return SKIPPAR;
    else if (strncmp(word, "BONDS", 5) == 0) return BONDPAR;
    else if (strncmp(word, "ANGLES", 6) == 0) return ANGLEPAR;
    else if (strncmp(word, "DIHEDRAL", 8) == 0) return TORSIONPAR;
    else if (strncmp(word, "IMPROPER", 8) == 0) return IMPRPAR;
    else if (strncmp(word, "CMAP", 4) == 0) return CMAPPAR;
    else if (strncmp(word, "NONBOND", 7) == 0 || strncmp(word, "cutnb", 5) == 0) return LJPAR;
    else if (strncmp(word, "HBOND", 5) == 0) return HBONDPAR;
    else if (strncmp(word, "END", 3) == 0) return ENDPAR;
    else
    {
        return CONTENTPAR;
    }
}

int skipLine(const char *line)
{
    if (line[0] == '\0' || line[0] == '\t' || line[0] == '\n' || line[0] == '!' || line[0] == '*')
    {
        return 1;
    }
    int count = 0;
    char * parse;
    char * copy = malloc(strlen(line) + 1);
    char * copy2 = malloc(strlen(line) + 1);
    strcpy(copy, line);
    strcpy(copy2, line);
    char match[] = "!";
    parse = strtok(copy, " ");
    while (parse != NULL)
    {
        int temp = strcmp(match, parse);
        if (temp == 0 && count == 0)
        {
            free(copy);
            free(copy2);
            return 1;
        }
        if (parse[0] == '!' && count == 0)
        {
            free(copy);
            free(copy2);
            return 1;
        }

        count++;
        parse = strtok(NULL, " ");
    }

    if (strcmp(copy2, "NBFIX\n") == 0)
    {
        free(copy);
        free(copy2);
        return 1;
    }

    free(copy);
    free(copy2);
    return 0;
}

int getParSize(const char *fileName, int* maxbondptr, int* maxangleptr, int* maxtorptr, int* maximpptr, int* maxljptr, int* maxcmapptr)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
    {
        printf("Failed to open file %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    *maxbondptr = 0;
    *maxangleptr = 0;
    *maxtorptr = 0;
    *maximpptr = 0;
    *maxljptr = 0;
    *maxcmapptr = 0;

    ParType secType = SKIPPAR;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        ParType lt = parType(line);
        if (lt == BONDPAR || lt == ANGLEPAR || lt == TORSIONPAR || lt == IMPRPAR ||
                lt == CMAPPAR || lt == LJPAR || lt == HBONDPAR)
        {
            secType = lt;
        }
        if (lt == CONTENTPAR)
        {
            if (secType == BONDPAR)
            {
                ++(*maxbondptr);
            }
            else if (secType == ANGLEPAR)
            {
                ++(*maxangleptr);
            }
            else if (secType == TORSIONPAR)
            {
                ++(*maxtorptr);
            }
            else if (secType == IMPRPAR)
            {
                ++(*maximpptr);
            }
            else if (secType == LJPAR)
            {
                ++(*maxljptr);
            }
            else if (secType == CMAPPAR)
            {
                *maxcmapptr = 6;
            }
        }
    }

    fclose(fp);
    if (line)
        free(line);

    return 0;

}

int parseCharmmPar(const char *fileName, CHARMM_PAR *charmmPars)
{

    int MAXBONDPAR = 0;
    int MAXANGLEPAR = 0;
    int MAXTORPAR = 0;
    int MAXIMPPAR = 0;
    int MAXLJPAR = 0;
    int MAXCMAPPAR = 0;

    getParSize(fileName, &MAXBONDPAR, &MAXANGLEPAR, &MAXTORPAR, &MAXIMPPAR, &MAXLJPAR, &MAXCMAPPAR);

    BOND_PARMS* bondStack = ddcMalloc(MAXBONDPAR * sizeof (BOND_PARMS));
    ANGLE_PARMS* angleStack = ddcMalloc(MAXANGLEPAR * sizeof (ANGLE_PARMS));
    TORSION_PARMS* torStack = ddcMalloc(MAXTORPAR * sizeof (TORSION_PARMS));
    IMPROPER_PARMS* imprStack = ddcMalloc(MAXIMPPAR * sizeof (IMPROPER_PARMS));
    LJCH_PARMS* ljStack = ddcMalloc(MAXLJPAR * sizeof (LJCH_PARMS));
    CMAP_PARMS* cmapStack = ddcMalloc(MAXCMAPPAR * sizeof (CMAP_PARMS));

    charmmPars->bondParmSize = MAXBONDPAR;
    charmmPars->bondParms = ddcMalloc(MAXBONDPAR * sizeof (BOND_PARMS*));

    charmmPars->angleParmSize = MAXANGLEPAR;
    charmmPars->angleParms = ddcMalloc(MAXANGLEPAR * sizeof (ANGLE_PARMS*));

    charmmPars->torParmSize = MAXTORPAR;
    charmmPars->torParms = ddcMalloc(MAXTORPAR * sizeof (TORSION_PARMS*));

    charmmPars->imprParmsize = MAXIMPPAR;
    charmmPars->imprParms = ddcMalloc(MAXIMPPAR * sizeof (IMPROPER_PARMS*));

    charmmPars->ljParmSize = MAXLJPAR;
    charmmPars->ljParms = ddcMalloc(MAXLJPAR * sizeof (LJCH_PARMS*));

    charmmPars->cmapParmSize = MAXCMAPPAR;
    charmmPars->cmapParms = ddcMalloc(MAXCMAPPAR * sizeof (CMAP_PARMS));

    int bondNum = 0;
    int angleNum = 0;
    int torNum = 0;
    int impNum = 0;
    int ljNum = 0;
    int cmapNum = 0;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fileName, "r");
    if (fp == NULL)
    {
        printf("Failed to open file %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    ParType secType = SKIPPAR;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        ParType lt = parType(line);
        if (lt == BONDPAR || lt == ANGLEPAR || lt == TORSIONPAR || lt == IMPRPAR ||
                lt == CMAPPAR || lt == LJPAR || lt == HBONDPAR)
        {
            secType = lt;
        }
        if (lt == CONTENTPAR)
        {
            // remove the comment and white space at the end of line
            char *noCLine = stripComment(line);
            char *newLine = rstrip(noCLine);

            if (secType == BONDPAR)
            {
                BOND_PARMS *bondParm = &bondStack[bondNum];
                if (processBondContent(newLine, bondParm) < 0)
                {
                    printf("Bond error: %s", line);
                }
                charmmPars->bondParms[bondNum] = bondParm;
                bondNum++;
            }
            else if (secType == ANGLEPAR)
            {
                ANGLE_PARMS *angleParm = &angleStack[angleNum];
                if (processAngleContent(newLine, angleParm) < 0)
                {
                    printf("Angle error: %s", line);
                }
                charmmPars->angleParms[angleNum] = angleParm;
                angleNum++;
            }
            else if (secType == TORSIONPAR)
            {
                TORSION_PARMS *torParm = &torStack[torNum];
                if (processTorContent(newLine, torParm) < 0)
                {
                    printf("Torsion error: %s", line);
                }
                charmmPars->torParms[torNum] = torParm;
                torNum++;
            }
            else if (secType == IMPRPAR)
            {
                IMPROPER_PARMS *imprParm = &imprStack[impNum];
                if (processImpContent(newLine, imprParm) < 0)
                {
                    printf("Angle error: %s", line);
                }
                charmmPars->imprParms[impNum] = imprParm;
                impNum++;
            }
            else if (secType == LJPAR)
            {
                LJCH_PARMS *ljParm = &ljStack[ljNum];
                if (processLJContent(newLine, ljParm) < 0)
                {
                    printf("Angle error: %s", line);
                }
                charmmPars->ljParms[ljNum] = ljParm;
                ljNum++;
            }
            else if (secType == CMAPPAR)
            {
                int gridSize = 24 * 24;
                int cmapID = cmapNum / gridSize;
                int gridIndex = cmapNum % gridSize;
                int nGridTerms[1];
                nGridTerms[0] = 0;
                CMAP_PARMS *cmapParm = &cmapStack[cmapID];
                if (processCMAPContent(newLine, cmapParm, gridIndex, nGridTerms) < 0)
                {
                    printf("CMAP error: %s", line);
                }
                charmmPars->cmapParms[cmapID] = cmapParm;
                cmapNum += nGridTerms[0];
            }
        }
    }
    fclose(fp);
    return 0;
}

int processBondContent(char* line, BOND_PARMS *bondParm)
{

    char *token;
    token = strtok(line, " ");
    int count = 0;

    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 0)
        {
            strncpy(bondParm->atmTypeI, token, 4);
            bondParm->atmTypeI[4] = '\0';
        }
        if (count == 1)
        {
            strncpy(bondParm->atmTypeJ, token, 4);
            bondParm->atmTypeJ[4] = '\0';
        }
        if (count == 2)
        {
            double kb;
            sscanf(token, "%lf", &kb);
            //bondParm->kb=units_convert(kb, "kcal*Angstrom^-2", NULL)/NA;
            bondParm->kb = units_convert(kb, "kcal*mol^-1*Angstrom^-2", NULL);
        }
        if (count == 3)
        {
            double b0;
            sscanf(token, "%lf", &b0);
            bondParm->b0 = units_convert(b0, "Angstrom", NULL);
            return 0;
        }

        token = strtok(NULL, " ");
        count++;
    }

    return -1;
}

int processAngleContent(char* line, ANGLE_PARMS *angleParm)
{

    char *token;
    token = strtok(line, " ");
    int count = 0;

    angleParm->kub = 0;
    angleParm->s0 = 0;
    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 0)
        {
            strncpy(angleParm->atmTypeI, token, 4);
            angleParm->atmTypeI[4] = '\0';
        }
        else if (count == 1)
        {
            strncpy(angleParm->atmTypeJ, token, 4);
            angleParm->atmTypeJ[4] = '\0';
        }
        else if (count == 2)
        {
            strncpy(angleParm->atmTypeK, token, 4);
            angleParm->atmTypeK[4] = '\0';
        }
        else if (count == 3)
        {
            double ktheta;
            sscanf(token, "%lf", &ktheta);
            //angleParm->ktheta=units_convert(ktheta, "kcal", NULL)/NA;
            angleParm->ktheta = units_convert(ktheta, "kcal*mol^-1", NULL);
        }
        else if (count == 4)
        {
            double theta0;
            sscanf(token, "%lf", &theta0);
            angleParm->theta0 = theta0 / RAD2DEG;
        }
        else if (count == 5)
        {
            double kub;
            sscanf(token, "%lf", &kub);
            //angleParm->kub=units_convert(kub, "kcal*Angstrom^-2", NULL)/NA;
            angleParm->kub = units_convert(kub, "kcal*mol^-1*Angstrom^-2", NULL);
        }
        else if (count == 6)
        {
            double s0;
            sscanf(token, "%lf", &s0);
            angleParm->s0 = units_convert(s0, "Angstrom", NULL);
        }

        token = strtok(NULL, " ");
        count++;
    }

    return 0;
}

int processTorContent(char* line, TORSION_PARMS *torParm)
{

    char *token;
    token = strtok(line, " ");
    int count = 0;

    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 0)
        {
            strncpy(torParm->atmTypeI, token, 4);
            torParm->atmTypeI[4] = '\0';
        }
        else if (count == 1)
        {
            strncpy(torParm->atmTypeJ, token, 4);
            torParm->atmTypeJ[4] = '\0';
        }
        else if (count == 2)
        {
            strncpy(torParm->atmTypeK, token, 4);
            torParm->atmTypeK[4] = '\0';
        }
        else if (count == 3)
        {
            strncpy(torParm->atmTypeL, token, 4);
            torParm->atmTypeL[4] = '\0';
        }
        else if (count == 4)
        {
            double kchi;
            sscanf(token, "%lf", &kchi);
            //torParm->kchi=units_convert(kchi, "kcal", NULL)/NA;
            torParm->kchi = units_convert(kchi, "kcal*mol^-1", NULL);
        }
        else if (count == 5)
        {
            sscanf(token, "%d", &torParm->n);
        }
        else if (count == 6)
        {
            double delta;
            sscanf(token, "%lf", &delta);
            torParm->delta = delta / RAD2DEG;
            return 0;
        }

        token = strtok(NULL, " ");
        count++;
    }

    return -1;
}

int processImpContent(char* line, IMPROPER_PARMS *imprParm)
{

    char *token;
    token = strtok(line, " ");
    int count = 0;

    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 0)
        {
            strncpy(imprParm->atmTypeI, token, 4);
            imprParm->atmTypeI[4] = '\0';
        }
        else if (count == 1)
        {
            strncpy(imprParm->atmTypeJ, token, 4);
            imprParm->atmTypeJ[4] = '\0';
        }
        else if (count == 2)
        {
            strncpy(imprParm->atmTypeK, token, 4);
            imprParm->atmTypeK[4] = '\0';
        }
        else if (count == 3)
        {
            strncpy(imprParm->atmTypeL, token, 4);
            imprParm->atmTypeL[4] = '\0';
        }
        else if (count == 4)
        {
            double kpsi;
            sscanf(token, "%lf", &kpsi);
            //imprParm->kpsi=units_convert(kpsi, "kcal", NULL)/NA;
            imprParm->kpsi = units_convert(kpsi, "kcal*mol^-1", NULL);
        }
        else if (count == 6)
        {
            double psi0;
            sscanf(token, "%lf", &psi0);
            imprParm->psi0 = psi0 / RAD2DEG;
            return 0;
        }

        token = strtok(NULL, " ");
        count++;
    }

    return -1;
}

int processLJContent(char* line, LJCH_PARMS *ljParm)
{

    char *token;
    token = strtok(line, " ");
    int count = 0;

    ljParm->eps14 = 0;
    ljParm->rmin14 = 0;
    /* walk through other tokens */
    while (token != NULL)
    {
        if (count == 0)
        {
            strncpy(ljParm->atmTypeI, token, 4);
            ljParm->atmTypeI[4] = '\0';
        }
        else if (count == 2)
        {
            double epsilon;
            sscanf(token, "%lf", &epsilon);
            //ljParm->epsilon=units_convert(epsilon, "kcal", NULL)/NA;
            ljParm->epsilon = units_convert(epsilon, "kcal*mol^-1", NULL);
        }
        else if (count == 3)
        {
            double rmin;
            sscanf(token, "%lf", &rmin);
            ljParm->rmin = units_convert(rmin, "Angstrom", NULL);
        }
        else if (count == 5)
        {
            double eps14;
            sscanf(token, "%lf", &eps14);
            //ljParm->eps14=units_convert(eps14, "kcal", NULL)/NA;
            ljParm->eps14 = units_convert(eps14, "kcal*mol^-1", NULL);
        }
        else if (count == 6)
        {
            double rmin14;
            sscanf(token, "%lf", &rmin14);
            ljParm->rmin14 = units_convert(rmin14, "Angstrom", NULL);
        }

        token = strtok(NULL, " ");
        count++;
    }

    return 0;
}

int processCMAPContent(char* line, CMAP_PARMS *cmapParm, int gridStart, int* nGridTerms)
{
    char *token;
    token = strtok(line, " \t\n");
    while (token != NULL)
    {
        if (line[0] != '!' && line[0] != 'C')
        {
            double value = strtod(token, NULL);
            int gridIndex = gridStart + nGridTerms[0];
            cmapParm->grid[gridIndex] = value;
            nGridTerms[0]++;
        }
        token = strtok(NULL, " \t\n");
    }
    return 0;
}

int freeCharmmPar(CHARMM_PAR *charmmPars)
{

    for (int i = 0; i < charmmPars->bondParmSize; i++)
    {
        free(charmmPars->bondParms[i]);
    }
    free(charmmPars->bondParms);

    for (int i = 0; i < charmmPars->angleParmSize; i++)
    {
        free(charmmPars->angleParms[i]);
    }
    free(charmmPars->angleParms);

    for (int i = 0; i < charmmPars->torParmSize; i++)
    {
        free(charmmPars->torParms[i]);
    }
    free(charmmPars->torParms);

    for (int i = 0; i < charmmPars->imprParmsize; i++)
    {
        free(charmmPars->imprParms[i]);
    }
    free(charmmPars->imprParms);

    for (int i = 0; i < charmmPars->ljParmSize; i++)
    {
        free(charmmPars->ljParms[i]);
    }
    free(charmmPars->ljParms);

    for (int i = 0; i < charmmPars->cmapParmSize; i++)
    {
        free(charmmPars->cmapParms[i]);
    }
    free(charmmPars->cmapParms);

    for (int i = 0; i < charmmPars->cmapParmSize; i++)
    {
        free(charmmPars->cmapParms[i]);
    }
    free(charmmPars->cmapParms);
    free(charmmPars);

    return 0;
}
