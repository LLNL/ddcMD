#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>

#include "bioCharmmTop.h"
#include "bioTransform.h"
#include "ddcMalloc.h"
//#include "cdbInterface.h"
//#include "mpiUtils.h"
#include "object.h"

int charmmTransformInit(BIOWEIGHTS *weights)
{
    weights->useFiles = 1;

    weights->weights = ddcMalloc(4 * sizeof (double));
    //Initialize  as 1
    for (int i = 0; i < 3; i++)
    {
        weights->weights[i] = 1.0;
    }
    //final is 0.0
    weights->weights[3] = 0.0;

    return 0;
}

int charmmTransformFree(BIOWEIGHTS *weights)
{
    ddcFree(weights->weights);
    return 0;
}

int readWeights(double *w1, double *w2)
{
    int err = 0;
    //  if(getRank(0) ==0) {
    //    char temp[512];
    char key[1000];
    struct stat buff;
    //getcwd(temp, 255);
    getcwd(key, 255);
    //    sprintf(temp+strlen(temp), "_%d_", loop);
    //    sprintf(key, "%s", temp);
    sprintf(key + strlen(key), "/weights.dat");

    printf("CHECK: readWeights file = %s\n", key);

    err = stat(key, &buff);
    if (err != -1)
    {
        FILE *file;
        char line[256];
        file = fopen(key, "r");
        if (file)
        {
            fgets(line, 255, file);
            printf("line = %s\n", line);
            //double temp = atof(line);
            //printf("temp + %f\n", temp);  
            *w1 = atof(line);
            //*w1 = temp;
            fgets(line, 255, file);
            printf("line = %s\n", line);
            //temp = atof(line);
            *w2 = atof(line);
            //*w2 =temp;
        }
        fclose(file);
    }
    printf("CHECK: w1 = %f w2 = %f\n", *w1, *w2);
    //} 
    //MPI_Bcast(&w1, 1, MPI_DOUBLE, 0, COMM_LOCAL);
    //MPI_Bcast(&w2, 1, MPI_DOUBLE, 0, COMM_LOCAL);
    return err;
}

//weights updated with cdb at the system loop level

int updateWeights(BIOWEIGHTS *weights)
{
    int err = 0;

    //  if(cdb->useCdb){
    //err = readTransformWeights(cdb, &weights->weights[1], &weights->weights[2], loop);

    //  } else {
    err = readWeights(&weights->weights[1], &weights->weights[2]);
    //  }
    //  printf("CHECK: read Weights w1= %f and w2 = %f\n", weights->weights[1], weights->weights[2]); 
    return err;
}

int getType(int atmI, int atmJ)
{
    if (atmI == 3 || atmJ == 3)
    {
        return 3;
    }
    else if (atmI == atmJ)
    {
        return atmI;
    }
    else if (atmI == 0)
    {
        return atmJ;
    }
    else if (atmJ == 0)
    {
        return atmI;
    }
    else
    {
        return 3;
    }
}

double get2Weights(BIOWEIGHTS *weights, TATOM_PARMS* atmI, TATOM_PARMS* atmJ)
{
    int speciesI = atmI->species;
    int speciesJ = atmJ->species;
    //printf("\tspeciesI = %d speciesJ = %d\n", speciesI, speciesJ);
    if (speciesI == speciesJ)
    {
        return weights->weights[speciesI];
    }
    else if ((speciesI == 0 && (speciesJ == 1 || speciesJ == 2)))
    {
        return weights->weights[speciesJ];
    }
    else if ((speciesJ == 0 && (speciesI == 1 || speciesI == 2)))
    {
        return weights->weights[speciesI];
    }
    else
    {
        return 0.0;
    }
}

double get3Weights(BIOWEIGHTS *weights, TATOM_PARMS* atmI, TATOM_PARMS* atmJ, TATOM_PARMS* atmK)
{
    int speciesI = atmI->species;
    int speciesJ = atmJ->species;
    int speciesK = atmK->species;

    int tempIJ = getType(speciesI, speciesJ);
    int typeIJK = getType(tempIJ, speciesK);

    //printf("CHECK get3Weights: (%d, %d, %d) -> index = %d weight = %f\n", speciesI, speciesJ, speciesK, typeIJK, weights->weights[typeIJK]);  

    return weights->weights[typeIJK];
}

double get4Weights(BIOWEIGHTS *weights, TATOM_PARMS* atmI, TATOM_PARMS* atmJ, TATOM_PARMS* atmK, TATOM_PARMS* atmL)
{
    int speciesI = atmI->species;
    int speciesJ = atmJ->species;
    int speciesK = atmK->species;
    int speciesL = atmL->species;

    int tempIJ = getType(speciesI, speciesJ);
    int tempKL = getType(speciesK, speciesL);
    int type = getType(tempIJ, tempKL);

    //printf("CHECK get4Weights: (%d, %d, %d, %d) -> index = %d weight = %f\n", speciesI, speciesJ, speciesK, speciesL, type, weights->weights[type]);  

    return weights->weights[type];
}

