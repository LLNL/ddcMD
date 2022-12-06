#include "bondedGPU.h"
#include "cudaUtils.h"
#include "units.h"
#include "ddcMath.h"
#include <thrust/pair.h>
#include "pairProcessGPU.h"
#include "bioCharmmParms.h"
#include <cmath>
#include <cub/cub.cuh>
#include <cub/device/device_radix_sort.cuh>
#include "bioGid.h"
#include "codata.h"
#include "gpu_allocator.hpp"
//#define CHECK(x) x  
#define CHECK(x)      // turning of logging
#define SHARED_BLIST 64  //shared mem perblock in nlist kernels
#define FLOAT_EPS    1e-08
#define NEAR_ZERO_ANGLE 0.017453292519943295 
#define NEAR_180_ANGLE 3.12413936106985 

#define PBC(x) x 
#define DEBUG_ENERGY(x) 

static int first;

/* Allocates memory for Residue structure (ie 'neighbor list') for
 *  all residues, and then sends structures from host to device's
 *  memory in a serialized manner.
 *  Function also converts data structures to struct of arrays for gpu
 */
void allocResiConManaged(SYSTEM *sys, CHARMMPOT_PARMS *parms, CHARMM_PARMS* charmmParms)
{
    printf("starting serialization of resi conn\n");
    STATE * state = sys->collection->state;
    int atomListSize = 0;
    int bondListSize = 0;
    int angleListSize = 0;
    int torsListSize = 0;

    gpu_allocator(parms->gpu_bondparms, 1);
    gpu_allocator(parms->gpu_angleparms, 1);
    gpu_allocator(parms->gpu_torsionparms, 1);

    CharmmBONDGPU_PARMS *gpu_parms = parms->gpu_bondparms;
    CharmmANGLEGPU_PARMS *gpu_angleparms = parms->gpu_angleparms;
    //CharmmTORSIONGPU_PARMS *gpu_torsparms = parms->gpu_torsionparms;

    //find buffer sizes needed for master bond/angle/torsion lists
    //by combining bond/angle/torsion lists of each residue
    CHECK(printf("get sizes\n");)
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        atomListSize += resiConn->atomListSize;
        bondListSize += resiConn->bondListSize;
        angleListSize += resiConn->angleListSize;
        torsListSize += resiConn->torsListSize;
    }

    //allocate memory for 'struct of array'
    CHECK(printf("begin managed allocation \n");)
    //atom info
    gpu_allocator(gpu_parms->numAtomsInResidue, charmmParms->resiConnSize);
    gpu_allocator(gpu_parms->numAtomsInResiduePrefixSum, charmmParms->resiConnSize);
    gpu_allocator(gpu_parms->rxs, sys->nion);
    gpu_allocator(gpu_parms->rys, sys->nion);
    gpu_allocator(gpu_parms->rzs, sys->nion);
    gpu_allocator(gpu_parms->resID, state->nion);
    gpu_allocator(gpu_parms->resIDs, state->nion);

    CHECK(printf("begin managed allocation for bonds\n");)
    //bond related buffers
    gpu_allocator(gpu_parms->bondListI, bondListSize);
    gpu_allocator(gpu_parms->bondListJ, bondListSize);
    gpu_allocator(gpu_parms->resiBondListStarts, charmmParms->resiConnSize);
    gpu_allocator(gpu_parms->resiBondListSizes, charmmParms->resiConnSize);
    gpu_allocator(gpu_parms->atmBondListStarts, atomListSize);
    gpu_allocator(gpu_parms->atmBondListSizes, atomListSize);
    gpu_allocator(gpu_parms->kb, bondListSize);
    gpu_allocator(gpu_parms->b0, bondListSize);

    CHECK(printf("begin managed allocation for angles\n");)
    //angle related buffers
    gpu_allocator(gpu_angleparms->angleListI, bondListSize);
    gpu_allocator(gpu_angleparms->angleListJ, bondListSize);
    gpu_allocator(gpu_angleparms->angleListK, bondListSize);
    gpu_allocator(gpu_angleparms->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_allocator(gpu_angleparms->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_allocator(gpu_angleparms->atmAngleListStarts, atomListSize);
    gpu_allocator(gpu_angleparms->atmAngleListSizes, atomListSize);
    gpu_allocator(gpu_angleparms->ktheta, angleListSize);
    gpu_allocator(gpu_angleparms->theta0, angleListSize);
    gpu_allocator(gpu_angleparms->kub, angleListSize);
    gpu_allocator(gpu_angleparms->s0, angleListSize);

    //load bond, angle, and tors lists into 1 unified master bond list to send to gpu
    int iat = 0; //total num atoms
    int ib = 0; //total bonds counter 
    int ian = 0; //total angles counter
    //int it = 0; //total torsions counter

    //for deach residue we, iterate over all bonds, angles, and torsions 
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        CHECK(printf("serializing residue %i\n", i);)
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        //printf("resi name %s \n",&(resiConn->resName) );
        BOND_CONN** bondList = resiConn->bondList;
        ANGLE_CONN** angleList = resiConn->angleList;
        //TORS_CONN** torsList = resiConn->torsList;

        //create gpu lists with atomlist info
        gpu_parms->numAtomsInResidue[i] = resiConn->atomListSize;
        //prefix sum atom list size
        gpu_parms->numAtomsInResiduePrefixSum[i] = iat;
        iat += resiConn->atomListSize;

        //create gpu lists with bond list info
        //1) where this residue's bond list begins in master bond list
        //2) the number bonds in this residue (where residue's bonds end in master bond list)
        gpu_parms->resiBondListStarts[i] = ib;
        gpu_parms->resiBondListSizes[i] = resiConn->bondListSize;
        //now for angles
        gpu_angleparms->resiAngleListStarts[i] = ian;
        gpu_angleparms->resiAngleListSizes[i] = resiConn->angleListSize;
        //torsions
        // gpu_torsparms->resiTorsionListStarts[i]=it;
        // gpu_torsparms->resiTorsionListSizes[i] = resiConn->torsListSize;

        printf("num bond %i\n", resiConn->bondListSize);
        printf("num angle %i\n", resiConn->angleListSize);
        printf("num tors %i\n", resiConn->torsListSize);

        //insert bonds for current residue into master bond list
        for (int j = 0; j < resiConn->bondListSize; j++)
        {
            gpu_parms->bondListI[ib] = bondList[j]->atmI;
            gpu_parms->bondListJ[ib] = bondList[j]->atmJ;
            gpu_parms->kb[ib] = bondList[j]->bondPtr->kb;
            gpu_parms->b0[ib] = bondList[j]->bondPtr->b0;
            ib++;
        }
        //now for angles
        for (int j = 0; j < resiConn->angleListSize; j++)
        {
            gpu_angleparms->angleListI[ian] = angleList[j]->atmI;
            gpu_angleparms->angleListJ[ian] = angleList[j]->atmJ;
            gpu_angleparms->angleListK[ian] = angleList[j]->atmK;
            gpu_angleparms->ktheta[ian] = angleList[j]->anglePtr->ktheta;
            gpu_angleparms->theta0[ian] = angleList[j]->anglePtr->theta0;
            gpu_angleparms->kub[ian] = angleList[j]->anglePtr->kub;
            gpu_angleparms->s0[ian] = angleList[j]->anglePtr->s0;
            ian++;
        }


        //create lists of bond, angle, and torsion list RANGES (per particle)
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            int atm_idx = gpu_parms->numAtomsInResiduePrefixSum[i] + j; //global 'atom type' index
            ATMRANGE *atmRange = resiConn->atmRanges[j];

            //bonds
            RANGE bondRange = atmRange->bondRange;
            if (bondRange.start == -1)
            {
                gpu_parms->atmBondListStarts[atm_idx] = -1;
                gpu_parms->atmBondListSizes[atm_idx] = 0;
            }
            else
            {
                int bond_idx = gpu_parms->resiBondListStarts[i] + bondRange.start;
                gpu_parms->atmBondListStarts[atm_idx] = bond_idx;
                gpu_parms->atmBondListSizes[atm_idx] = bondRange.end - bondRange.start;
            }

            //angles
            RANGE angleRange = atmRange->angleRange;
            if (angleRange.start == -1)
            {
                gpu_angleparms->atmAngleListStarts[atm_idx] = -1;
                gpu_angleparms->atmAngleListSizes[atm_idx] = 0;
            }
            else
            {
                int angle_idx = gpu_angleparms->resiAngleListStarts[i] + angleRange.start;
                gpu_angleparms->atmAngleListStarts[atm_idx] = angle_idx;
                gpu_angleparms->atmAngleListSizes[atm_idx] = angleRange.end - angleRange.start;
            }
        }
        /*
                   //load angle list
                   ANGLE_CONN** angleList = resiConn->angleList;
                   for (int j =0 ; j < resiConn->angleListSize; j++)
                   {
                           angleListI[ia] = angleList[j]->atmI;
                           angleListJ[ia] = angleList[j]->atmJ;
                           angleListK[ia] = angleList[j]->atmK;
                           ia++;
                   }

                   //load torsion list
                   TORS_CONN** torsList = resiConn->torsList;
                   for (int j =0 ; j < resiConn->torsListSize; j++)
                   {
                           torsListI[it] = torsList[j]->atmI;
                           torsListJ[it] = torsList[j]->atmJ;
                           torsListK[it] = torsList[j]->atmK;
                           torsListL[it] = torsList[j]->atmL;
                           it++;
                   }
         */
    }

    /*
    for (int i =0; i <atomListSize ; i++)
    {
       printf(" %i bond start %i size %i \n",i, gpu_parms->atmBondListStarts[i], gpu_parms->atmBondListSizes[i]);
       printf(" %i angle start %i size %i \n",i, gpu_angleparms->atmAngleListStarts[i], gpu_angleparms->atmAngleListSizes[i]);
    }

    for (int i =0; i< bondListSize; i++){
       printf("i %i\n", i);
       printf("b0 %f \n", gpu_parms->b0[i]);
       printf("kb %f \n", gpu_parms->kb[i]);
    }

    for (int i =0; i< angleListSize; i++){
       printf("i %i\n", i);
       printf("kub %f \n", gpu_angleparms->kub[i]);
       printf("s0 %f \n", gpu_angleparms->s0[i]);
    }
     */
    //assign residue type ids to each particle
    for (int i = 0; i < state->nion; i++)
    {
        char *name = state->species[i]->name;
        RESI_CONN * resi_conn = findResiConnNew(charmmParms, name);
        //printf("particle %i %s %i \n", i, name, resi_conn->resID-1);
        gpu_parms->resID[i] = resi_conn->resID - 1;
    }

    printf("bl %i angle %i tors %i\n", bondListSize, angleListSize, torsListSize);
}

/*
To the poor, tortured soul having to read this function:
Yes, I tried using unified memory. It resulted in weird heap overflow errors
Hence the reason for this monstrosity. I am truly sorry for any suffering this causes
This code does not reflect my standards for code quality. 
 */
void allocResiCon(SYSTEM *sys, CHARMMPOT_PARMS *parms) 
{
    GPUNLIST *gnlist=sys->collection->gnlist;
    STATE * state = sys->collection->state;
    first = 1;
    CHECK(printf("begin serialized allocations \n");)
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    
    int atomListSize = 0;
    int bondListSize = 0;
    int angleListSize = 0;
    int cosangleListSize = 0;
    int rebangleListSize = 0;
    int torsionListSize = 0;
    int improperListSize = 0;
    int bpairListSize = 0;

    printf("nion sys %i state %i\n", sys->nion, state->nion);
    //alloc temp host struct for serializing
    CharmmBONDGPU_PARMS* bondparms_h = (CharmmBONDGPU_PARMS*) malloc(sizeof (CharmmBONDGPU_PARMS));
    CharmmANGLEGPU_PARMS* angleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    CharmmANGLEGPU_PARMS* cosangleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    CharmmANGLEGPU_PARMS* rebangleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    CharmmTORSIONGPU_PARMS* torsionparms_h = (CharmmTORSIONGPU_PARMS*) malloc(sizeof (CharmmTORSIONGPU_PARMS));
    CharmmIMPROPERGPU_PARMS* improperparms_h = (CharmmIMPROPERGPU_PARMS*) malloc(sizeof (CharmmIMPROPERGPU_PARMS));
    CharmmBPAIRGPU_PARMS* bpairparms_h = (CharmmBPAIRGPU_PARMS*) malloc(sizeof (CharmmBPAIRGPU_PARMS));

    parms->bondparms_h = bondparms_h;
    parms->angleparms_h = angleparms_h;
    parms->cosangleparms_h = cosangleparms_h;
    parms->rebangleparms_h = rebangleparms_h;
    parms->torsionparms_h = torsionparms_h;
    parms->improperparms_h = improperparms_h;
    parms->bpairparms_h = bpairparms_h;


    //allocate host struct with pointers to gpu mem
    parms->gpu_bondparms_h = (CharmmBONDGPU_PARMS*) malloc(sizeof (CharmmBONDGPU_PARMS));
    parms->gpu_angleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    parms->gpu_cosangleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    parms->gpu_rebangleparms_h = (CharmmANGLEGPU_PARMS*) malloc(sizeof (CharmmANGLEGPU_PARMS));
    parms->gpu_torsionparms_h = (CharmmTORSIONGPU_PARMS*) malloc(sizeof (CharmmTORSIONGPU_PARMS));
    parms->gpu_improperparms_h = (CharmmIMPROPERGPU_PARMS*) malloc(sizeof (CharmmIMPROPERGPU_PARMS));
    parms->gpu_bpairparms_h = (CharmmBPAIRGPU_PARMS*) malloc(sizeof (CharmmBPAIRGPU_PARMS));

    //allocate device structs with pointers to gpu mem
    gpu_allocator(parms->gpu_bondparms, 1);
    gpu_allocator(parms->gpu_angleparms, 1);
    gpu_allocator(parms->gpu_cosangleparms, 1);
    gpu_allocator(parms->gpu_rebangleparms, 1);
    gpu_allocator(parms->gpu_torsionparms, 1);
    gpu_allocator(parms->gpu_improperparms, 1);
    gpu_allocator(parms->gpu_bpairparms, 1);

    //find buffer sizes needed for master bond/angle/torsion lists
    //by combining bond/angle/torsion lists of each residue
    CHECK(printf("get sizes\n");)
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        atomListSize += resiConn->atomListSize;
        bondListSize += resiConn->bondListSize;
        angleListSize += resiConn->angleListSize;
        cosangleListSize += resiConn->cosangleListSize;
        rebangleListSize += resiConn->rebangleListSize;
        torsionListSize += resiConn->torsListSize;
        improperListSize += resiConn->imprListSize;
        bpairListSize += resiConn->bpairListSize;
    }

    CHECK(printf("commit sizes\n");)
    //commit counts to relevant structs
    bondparms_h->atomListSize = atomListSize;
    bondparms_h->bondListSize = bondListSize;
    angleparms_h->angleListSize = angleListSize;
    cosangleparms_h->angleListSize = cosangleListSize;
    rebangleparms_h->angleListSize = rebangleListSize;
    torsionparms_h->torsionListSize = torsionListSize;
    improperparms_h->improperListSize = improperListSize;
    bpairparms_h->bpairListSize = bpairListSize;

    parms->gpu_bondparms_h->atomListSize = atomListSize;
    parms->gpu_bondparms_h->bondListSize = bondListSize;
    parms->gpu_angleparms_h->angleListSize = angleListSize;
    parms->gpu_cosangleparms_h->angleListSize = cosangleListSize;
    parms->gpu_rebangleparms_h->angleListSize = rebangleListSize;
    parms->gpu_torsionparms_h->torsionListSize = torsionListSize;
    parms->gpu_improperparms_h->improperListSize = improperListSize;
    parms->gpu_bpairparms_h->bpairListSize = bpairListSize;


    CHECK(printf("alloc particle buffers sorted by gid\n");)
    //allocate postiion info. might want to move this for abstraction sake
    //only need to this on gpu because this doesn't need to get serialized
    gpu_allocator(parms->gpu_bondparms_h->rxs, sys->nion);
    gpu_allocator(parms->gpu_bondparms_h->rys, sys->nion);
    gpu_allocator(parms->gpu_bondparms_h->rzs, sys->nion);

    //allocate memory for 'struct of array'
    //atom info
    CHECK(printf("alloc atom info\n");)
    gpu_allocator(parms->gpu_bondparms_h->numAtomsInResidue, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_bondparms_h->numAtomsInResiduePrefixSum, charmmParms->resiConnSize);
    bondparms_h->numAtomsInResidue = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    bondparms_h->numAtomsInResiduePrefixSum = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);

    //box
    bondparms_h->h1 = gnlist->hmat_g;

    //residue info
    CHECK(printf("alloc resi info\n");)
    gpu_allocator(parms->gpu_bondparms_h->resID, sys->nion);
    gpu_allocator(parms->gpu_bondparms_h->resIDs, sys->nion);
    bondparms_h->resID = (int*) malloc(sizeof (int)*sys->nion);
    bondparms_h->resIDs = (int*) malloc(sizeof (int)*sys->nion);


    CHECK(printf("allocation for bonds gpu\n");)
    gpu_allocator(parms->gpu_bondparms_h->bondListI, bondListSize);
    gpu_allocator(parms->gpu_bondparms_h->bondListJ, bondListSize);
    gpu_allocator(parms->gpu_bondparms_h->resiBondListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_bondparms_h->resiBondListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_bondparms_h->atmBondListStarts, atomListSize);
    gpu_allocator(parms->gpu_bondparms_h->atmBondListSizes, atomListSize);
    gpu_allocator(parms->gpu_bondparms_h->kb, bondListSize);
    gpu_allocator(parms->gpu_bondparms_h->b0, bondListSize);


    CHECK(printf("allocation for bonds host\n");)
    bondparms_h->bondListI = (int*) malloc(sizeof (int)*bondListSize);
    bondparms_h->bondListJ = (int*) malloc(sizeof (int)*bondListSize);
    bondparms_h->resiBondListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    bondparms_h->resiBondListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    bondparms_h->atmBondListStarts = (int*) malloc(sizeof (int)*atomListSize);
    bondparms_h->atmBondListSizes = (int*) malloc(sizeof (int)*atomListSize);
    bondparms_h->kb = (double*) malloc(sizeof (double)*bondListSize);
    bondparms_h->b0 = (double*) malloc(sizeof (double)*bondListSize);


    CHECK(printf("allocation for angles gpu\n");)
    gpu_allocator(parms->gpu_angleparms_h->angleListI, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->angleListJ, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->angleListK, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_angleparms_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_angleparms_h->atmAngleListStarts, atomListSize);
    gpu_allocator(parms->gpu_angleparms_h->atmAngleListSizes, atomListSize);
    gpu_allocator(parms->gpu_angleparms_h->ktheta, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->theta0, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->kub, angleListSize);
    gpu_allocator(parms->gpu_angleparms_h->s0, angleListSize);


    CHECK(printf("allocation for angles host\n");)
    angleparms_h->angleListI = (int*) malloc(sizeof (int)*angleListSize);
    angleparms_h->angleListJ = (int*) malloc(sizeof (int)*angleListSize);
    angleparms_h->angleListK = (int*) malloc(sizeof (int)*angleListSize);
    angleparms_h->resiAngleListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    angleparms_h->resiAngleListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    angleparms_h->atmAngleListStarts = (int*) malloc(sizeof (int)*atomListSize);
    angleparms_h->atmAngleListSizes = (int*) malloc(sizeof (int)*atomListSize);
    angleparms_h->ktheta = (double*) malloc(sizeof (double)*angleListSize);
    angleparms_h->theta0 = (double*) malloc(sizeof (double)*angleListSize);
    angleparms_h->kub = (double*) malloc(sizeof (double)*angleListSize);
    angleparms_h->s0 = (double*) malloc(sizeof (double)*angleListSize);


    CHECK(printf("allocation for cosangles gpu\n");)
    gpu_allocator(parms->gpu_cosangleparms_h->angleListI, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->angleListJ, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->angleListK, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_cosangleparms_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_cosangleparms_h->atmAngleListStarts, atomListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->atmAngleListSizes, atomListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->ktheta, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->theta0, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->kub, cosangleListSize);
    gpu_allocator(parms->gpu_cosangleparms_h->s0, cosangleListSize);


    CHECK(printf("allocation for cos angles host\n");)
    cosangleparms_h->angleListI = (int*) malloc(sizeof (int)*cosangleListSize);
    cosangleparms_h->angleListJ = (int*) malloc(sizeof (int)*cosangleListSize);
    cosangleparms_h->angleListK = (int*) malloc(sizeof (int)*cosangleListSize);
    cosangleparms_h->resiAngleListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    cosangleparms_h->resiAngleListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    cosangleparms_h->atmAngleListStarts = (int*) malloc(sizeof (int)*atomListSize);
    cosangleparms_h->atmAngleListSizes = (int*) malloc(sizeof (int)*atomListSize);
    cosangleparms_h->ktheta = (double*) malloc(sizeof (double)*cosangleListSize);
    cosangleparms_h->theta0 = (double*) malloc(sizeof (double)*cosangleListSize);
    cosangleparms_h->kub = (double*) malloc(sizeof (double)*cosangleListSize);
    cosangleparms_h->s0 = (double*) malloc(sizeof (double)*cosangleListSize);

    CHECK(printf("allocation for restricted bending angles gpu\n");)
    gpu_allocator(parms->gpu_rebangleparms_h->angleListI, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->angleListJ, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->angleListK, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_rebangleparms_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_rebangleparms_h->atmAngleListStarts, atomListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->atmAngleListSizes, atomListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->ktheta, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->theta0, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->kub, rebangleListSize);
    gpu_allocator(parms->gpu_rebangleparms_h->s0, rebangleListSize);


    CHECK(printf("allocation for restricted bending angles host\n");)
    rebangleparms_h->angleListI = (int*) malloc(sizeof (int)*rebangleListSize);
    rebangleparms_h->angleListJ = (int*) malloc(sizeof (int)*rebangleListSize);
    rebangleparms_h->angleListK = (int*) malloc(sizeof (int)*rebangleListSize);
    rebangleparms_h->resiAngleListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    rebangleparms_h->resiAngleListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    rebangleparms_h->atmAngleListStarts = (int*) malloc(sizeof (int)*atomListSize);
    rebangleparms_h->atmAngleListSizes = (int*) malloc(sizeof (int)*atomListSize);
    rebangleparms_h->ktheta = (double*) malloc(sizeof (double)*rebangleListSize);
    rebangleparms_h->theta0 = (double*) malloc(sizeof (double)*rebangleListSize);
    rebangleparms_h->kub = (double*) malloc(sizeof (double)*rebangleListSize);
    rebangleparms_h->s0 = (double*) malloc(sizeof (double)*rebangleListSize);    

    CHECK(printf("allocation for torsion gpu\n");)
    gpu_allocator(parms->gpu_torsionparms_h->torsionListI, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->torsionListJ, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->torsionListK, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->torsionListL, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->resiTorsionListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_torsionparms_h->resiTorsionListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_torsionparms_h->atmTorsionListStarts, atomListSize);
    gpu_allocator(parms->gpu_torsionparms_h->atmTorsionListSizes, atomListSize);
    gpu_allocator(parms->gpu_torsionparms_h->kchi, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->delta, torsionListSize);
    gpu_allocator(parms->gpu_torsionparms_h->n, torsionListSize);


    CHECK(printf("allocation for torsions host\n");)
    torsionparms_h->torsionListI = (int*) malloc(sizeof (int)*torsionListSize);
    torsionparms_h->torsionListJ = (int*) malloc(sizeof (int)*torsionListSize);
    torsionparms_h->torsionListK = (int*) malloc(sizeof (int)*torsionListSize);
    torsionparms_h->torsionListL = (int*) malloc(sizeof (int)*torsionListSize);
    torsionparms_h->resiTorsionListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    torsionparms_h->resiTorsionListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    torsionparms_h->atmTorsionListStarts = (int*) malloc(sizeof (int)*atomListSize);
    torsionparms_h->atmTorsionListSizes = (int*) malloc(sizeof (int)*atomListSize);
    torsionparms_h->kchi = (double*) malloc(sizeof (double)*torsionListSize);
    torsionparms_h->delta = (double*) malloc(sizeof (double)*torsionListSize);
    torsionparms_h->n = (int*) malloc(sizeof (int)*torsionListSize);


    CHECK(printf("allocation for improper gpu\n");)
    gpu_allocator(parms->gpu_improperparms_h->improperListI, improperListSize);
    gpu_allocator(parms->gpu_improperparms_h->improperListJ, improperListSize);
    gpu_allocator(parms->gpu_improperparms_h->improperListK, improperListSize);
    gpu_allocator(parms->gpu_improperparms_h->improperListL, improperListSize);
    gpu_allocator(parms->gpu_improperparms_h->resiImproperListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_improperparms_h->resiImproperListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_improperparms_h->atmImproperListStarts, atomListSize);
    gpu_allocator(parms->gpu_improperparms_h->atmImproperListSizes, atomListSize);
    gpu_allocator(parms->gpu_improperparms_h->kpsi, improperListSize);
    gpu_allocator(parms->gpu_improperparms_h->psi0, improperListSize);


    CHECK(printf("allocation for impropers host\n");)
    improperparms_h->improperListI = (int*) malloc(sizeof (int)*improperListSize);
    improperparms_h->improperListJ = (int*) malloc(sizeof (int)*improperListSize);
    improperparms_h->improperListK = (int*) malloc(sizeof (int)*improperListSize);
    improperparms_h->improperListL = (int*) malloc(sizeof (int)*improperListSize);
    improperparms_h->resiImproperListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    improperparms_h->resiImproperListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    improperparms_h->atmImproperListStarts = (int*) malloc(sizeof (int)*atomListSize);
    improperparms_h->atmImproperListSizes = (int*) malloc(sizeof (int)*atomListSize);
    improperparms_h->kpsi = (double*) malloc(sizeof (double)*improperListSize);
    improperparms_h->psi0 = (double*) malloc(sizeof (double)*improperListSize);



    CHECK(printf("allocation for bpair gpu\n");)
    gpu_allocator(parms->gpu_bpairparms_h->bpairListI, bpairListSize);
    gpu_allocator(parms->gpu_bpairparms_h->bpairListJ, bpairListSize);
    gpu_allocator(parms->gpu_bpairparms_h->resiBpairListStarts, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_bpairparms_h->resiBpairListSizes, charmmParms->resiConnSize);
    gpu_allocator(parms->gpu_bpairparms_h->atmBpairListStarts, atomListSize);
    gpu_allocator(parms->gpu_bpairparms_h->atmBpairListSizes, atomListSize);
    gpu_allocator(parms->gpu_bpairparms_h->shift, bpairListSize);
    gpu_allocator(parms->gpu_bpairparms_h->sigma, bpairListSize);
    gpu_allocator(parms->gpu_bpairparms_h->eps, bpairListSize);
    gpu_allocator(parms->gpu_bpairparms_h->q, sys->nion);


    CHECK(printf("allocation for bpair host\n");)
    bpairparms_h->bpairListI = (int*) malloc(sizeof (int)*bpairListSize);
    bpairparms_h->bpairListJ = (int*) malloc(sizeof (int)*bpairListSize);
    bpairparms_h->resiBpairListStarts = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    bpairparms_h->resiBpairListSizes = (int*) malloc(sizeof (int)*charmmParms->resiConnSize);
    bpairparms_h->atmBpairListStarts = (int*) malloc(sizeof (int)*atomListSize);
    bpairparms_h->atmBpairListSizes = (int*) malloc(sizeof (int)*atomListSize);
    bpairparms_h->shift = (double*) malloc(sizeof (double)*bpairListSize);
    bpairparms_h->sigma = (double*) malloc(sizeof (double)*bpairListSize);
    bpairparms_h->eps = (double*) malloc(sizeof (double)*bpairListSize);

    //begin serialization
    //load bond, angle, and tors lists into 1 unified master bond list to send to gpu
    int iat = 0; //total num atoms
    int ib = 0; //total bonds counter 
    int ian = 0; //total angles counter
    int ican = 0; //total cosine angles counter
    int irban = 0; // total restrict bending angle counter
    int it = 0; //total torsions counter
    int iim = 0; //total impropers counter
    int ibp = 0; //total bpair counter

    //for deach residue we, iterate over all bonds, angles, and torsions 
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        CHECK(printf("serializing residue %i\n", i);)
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        //printf("resi name %s \n",&(resiConn->resName) );
        BOND_CONN** bondList = resiConn->bondList;
        ANGLE_CONN** angleList = resiConn->angleList;
        ANGLE_CONN** cosangleList = resiConn->cosangleList;
        ANGLE_CONN** rebangleList = resiConn->rebangleList;
        TORS_CONN** torsList = resiConn->torsList;
        IMPR_CONN** imprList = resiConn->imprList;
        BPAIR_CONN** bpairList = resiConn->bpairList;

        //create gpu lists with atomlist info
        bondparms_h->numAtomsInResidue[i] = resiConn->atomListSize;
        //prefix sum atom list size
        bondparms_h->numAtomsInResiduePrefixSum[i] = iat;
        iat += resiConn->atomListSize;

        //create gpu lists with bond list info
        //1) where this residue's bond list begins in master bond list
        //2) the number bonds in this residue (where residue's bonds end in master bond list)
        bondparms_h->resiBondListStarts[i] = ib;
        bondparms_h->resiBondListSizes[i] = resiConn->bondListSize;
        //now for angles
        angleparms_h->resiAngleListStarts[i] = ian;
        angleparms_h->resiAngleListSizes[i] = resiConn->angleListSize;

        //cosine angles
        cosangleparms_h->resiAngleListStarts[i] = ican;
        cosangleparms_h->resiAngleListSizes[i] = resiConn->cosangleListSize;

        //restricted bending angles
        rebangleparms_h->resiAngleListStarts[i] = irban;
        rebangleparms_h->resiAngleListSizes[i] = resiConn->rebangleListSize;
        
        //torsion
        torsionparms_h->resiTorsionListStarts[i] = it;
        torsionparms_h->resiTorsionListSizes[i] = resiConn->torsListSize;

        //improper
        improperparms_h->resiImproperListStarts[i] = iim;
        improperparms_h->resiImproperListSizes[i] = resiConn->imprListSize;

        //bpair
        bpairparms_h->resiBpairListStarts[i] = ibp;
        bpairparms_h->resiBpairListSizes[i] = resiConn->bpairListSize;
        /*
        printf("num bond %i\n", resiConn->bondListSize);
        printf("num angle %i\n", resiConn->angleListSize);
        printf("num cos angle %i\n", resiConn->cosangleListSize);
        printf("num tors %i\n", resiConn->torsListSize);
        printf("num impr %i\n", resiConn->imprListSize);
        printf("num bpair %i\n", resiConn->bpairListSize);
         */
        //insert bonds for current residue into master bond list
        for (int j = 0; j < resiConn->bondListSize; j++)
        {
            bondparms_h->bondListI[ib] = bondList[j]->atmI;
            bondparms_h->bondListJ[ib] = bondList[j]->atmJ;
            bondparms_h->kb[ib] = bondList[j]->bondPtr->kb;
            bondparms_h->b0[ib] = bondList[j]->bondPtr->b0;
            ib++;
        }

        //now for angles
        for (int j = 0; j < resiConn->angleListSize; j++)
        {
            angleparms_h->angleListI[ian] = angleList[j]->atmI;
            angleparms_h->angleListJ[ian] = angleList[j]->atmJ;
            angleparms_h->angleListK[ian] = angleList[j]->atmK;
            angleparms_h->ktheta[ian] = angleList[j]->anglePtr->ktheta;
            angleparms_h->theta0[ian] = angleList[j]->anglePtr->theta0;
            angleparms_h->kub[ian] = angleList[j]->anglePtr->kub;
            angleparms_h->s0[ian] = angleList[j]->anglePtr->s0;
            ian++;
        }

        //now for cos angles
        for (int j = 0; j < resiConn->cosangleListSize; j++)
        {
            cosangleparms_h->angleListI[ican] = cosangleList[j]->atmI;
            cosangleparms_h->angleListJ[ican] = cosangleList[j]->atmJ;
            cosangleparms_h->angleListK[ican] = cosangleList[j]->atmK;
            cosangleparms_h->ktheta[ican] = cosangleList[j]->anglePtr->ktheta;
            cosangleparms_h->theta0[ican] = cosangleList[j]->anglePtr->theta0;
            cosangleparms_h->kub[ican] = cosangleList[j]->anglePtr->kub;
            cosangleparms_h->s0[ican] = cosangleList[j]->anglePtr->s0;
            ican++;
        }

        //now for restricted bending angles
        for (int j = 0; j < resiConn->rebangleListSize; j++)
        {
            rebangleparms_h->angleListI[irban] = rebangleList[j]->atmI;
            rebangleparms_h->angleListJ[irban] = rebangleList[j]->atmJ;
            rebangleparms_h->angleListK[irban] = rebangleList[j]->atmK;
            rebangleparms_h->ktheta[irban] = rebangleList[j]->anglePtr->ktheta;
            rebangleparms_h->theta0[irban] = rebangleList[j]->anglePtr->theta0;
            rebangleparms_h->kub[irban] = rebangleList[j]->anglePtr->kub;
            rebangleparms_h->s0[irban] = rebangleList[j]->anglePtr->s0;
            irban++;
        }        
        
        //now for torsions
        for (int j = 0; j < resiConn->torsListSize; j++)
        {
            torsionparms_h->torsionListI[it] = torsList[j]->atmI;
            torsionparms_h->torsionListJ[it] = torsList[j]->atmJ;
            torsionparms_h->torsionListK[it] = torsList[j]->atmK;
            torsionparms_h->torsionListL[it] = torsList[j]->atmL;
            torsionparms_h->kchi[it] = torsList[j]->torsPtr->kchi;
            torsionparms_h->delta[it] = torsList[j]->torsPtr->delta;
            torsionparms_h->n[it] = torsList[j]->torsPtr->n;
            it++;
        }

        //now for impropers
        for (int j = 0; j < resiConn->imprListSize; j++)
        {
            improperparms_h->improperListI[iim] = imprList[j]->atmI;
            improperparms_h->improperListJ[iim] = imprList[j]->atmJ;
            improperparms_h->improperListK[iim] = imprList[j]->atmK;
            improperparms_h->improperListL[iim] = imprList[j]->atmL;
            improperparms_h->kpsi[iim] = imprList[j]->imprPtr->kpsi;
            improperparms_h->psi0[iim] = imprList[j]->imprPtr->psi0;
            iim++;
        }


        //bpair
        for (int j = 0; j < resiConn->bpairListSize; j++)
        {
            bpairparms_h->bpairListI[ibp] = bpairList[j]->atmI;
            bpairparms_h->bpairListJ[ibp] = bpairList[j]->atmJ;
            bpairparms_h->sigma[ibp] = bpairList[j]->sigma;
            bpairparms_h->eps[ibp] = bpairList[j]->eps;
            bpairparms_h->shift[ibp] = bpairList[j]->shift;
            ibp++;
        }

        //create lists of bond, angle, and torsion list RANGES (per particle)
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            int atm_idx = bondparms_h->numAtomsInResiduePrefixSum[i] + j; //global 'atom type' index
            ATMRANGE *atmRange = resiConn->atmRanges[j];

            //bonds
            RANGE bondRange = atmRange->bondRange;
            if (bondRange.start == -1)
            {
                bondparms_h->atmBondListStarts[atm_idx] = -1;
                bondparms_h->atmBondListSizes[atm_idx] = 0;
            }
            else
            {
                int bond_idx = bondparms_h->resiBondListStarts[i] + bondRange.start;
                bondparms_h->atmBondListStarts[atm_idx] = bond_idx;
                bondparms_h->atmBondListSizes[atm_idx] = bondRange.end - bondRange.start;
            }

            //angles
            RANGE angleRange = atmRange->angleRange;
            if (angleRange.start == -1)
            {
                angleparms_h->atmAngleListStarts[atm_idx] = -1;
                angleparms_h->atmAngleListSizes[atm_idx] = 0;
            }
            else
            {
                int angle_idx = angleparms_h->resiAngleListStarts[i] + angleRange.start;
                angleparms_h->atmAngleListStarts[atm_idx] = angle_idx;
                angleparms_h->atmAngleListSizes[atm_idx] = angleRange.end - angleRange.start;
            }

            //cos angles
            RANGE cosangleRange = atmRange->cosangleRange;
            if (cosangleRange.start == -1)
            {
                cosangleparms_h->atmAngleListStarts[atm_idx] = -1;
                cosangleparms_h->atmAngleListSizes[atm_idx] = 0;
            }
            else
            {
                int angle_idx = cosangleparms_h->resiAngleListStarts[i] + cosangleRange.start;
                cosangleparms_h->atmAngleListStarts[atm_idx] = angle_idx;
                cosangleparms_h->atmAngleListSizes[atm_idx] = cosangleRange.end - cosangleRange.start;
            }

            //restricted bending angles
            RANGE rebangleRange = atmRange->rebangleRange;
            if (rebangleRange.start == -1)
            {
                rebangleparms_h->atmAngleListStarts[atm_idx] = -1;
                rebangleparms_h->atmAngleListSizes[atm_idx] = 0;
            }
            else
            {
                int angle_idx = rebangleparms_h->resiAngleListStarts[i] + rebangleRange.start;
                rebangleparms_h->atmAngleListStarts[atm_idx] = angle_idx;
                rebangleparms_h->atmAngleListSizes[atm_idx] = rebangleRange.end - rebangleRange.start;
            }            
            
            //tors
            RANGE torsionRange = atmRange->torsRange;
            if (torsionRange.start == -1)
            {
                torsionparms_h->atmTorsionListStarts[atm_idx] = -1;
                torsionparms_h->atmTorsionListSizes[atm_idx] = 0;
            }
            else
            {
                int torsion_idx = torsionparms_h->resiTorsionListStarts[i] + torsionRange.start;
                torsionparms_h->atmTorsionListStarts[atm_idx] = torsion_idx;
                torsionparms_h->atmTorsionListSizes[atm_idx] = torsionRange.end - torsionRange.start;
            }

            //improper
            RANGE improperRange = atmRange->imprRange;
            if (improperRange.start == -1)
            {
                improperparms_h->atmImproperListStarts[atm_idx] = -1;
                improperparms_h->atmImproperListSizes[atm_idx] = 0;
            }
            else
            {
                int improper_idx = improperparms_h->resiImproperListStarts[i] + improperRange.start;
                improperparms_h->atmImproperListStarts[atm_idx] = improper_idx;
                improperparms_h->atmImproperListSizes[atm_idx] = improperRange.end - improperRange.start;
            }


            //bpairs
            RANGE bpairRange = atmRange->bpairRange;
            if (bpairRange.start == -1)
            {
                bpairparms_h->atmBpairListStarts[atm_idx] = -1;
                bpairparms_h->atmBpairListSizes[atm_idx] = 0;
            }
            else
            {
                int bpair_idx = bpairparms_h->resiBpairListStarts[i] + bpairRange.start;
                bpairparms_h->atmBpairListStarts[atm_idx] = bpair_idx;
                bpairparms_h->atmBpairListSizes[atm_idx] = bpairRange.end - bpairRange.start;
            }



        }
    }

    for (int i = 0; i < atomListSize; i++)
    {
        //printf(" %i tors start %i size %i \n",i, bondparms_h->atmBondListStarts[i], bondparms_h->atmBondListSizes[i]);
        //printf(" %i angle start %i size %i \n",i, angleparms_h->atmAngleListStarts[i], angleparms_h->atmAngleListSizes[i]);
        //printf(" %i tors start %i size %i \n",i, torsionparms_h->atmTorsionListStarts[i], torsionparms_h->atmTorsionListSizes[i]);
    }
    /*
         printf("bonds\n");
    for (int i =0; i< bondListSize; i++){
       printf("i %i b0 %f kb %f \n", i,bondparms_h->b0[i],bondparms_h->kb[i]);
    }

         printf("angles\n");
    for (int i =0; i< angleListSize; i++){
       printf("i %i\n", i);
       printf("kub %f \n", angleparms_h->kub[i]);
       printf("s0 %f \n", angleparms_h->s0[i]);
    }

         printf("cos angles\n");
    for (int i =0; i< cosangleListSize; i++){
       printf("i %i\n", i);
       printf("kub %f \n", cosangleparms_h->kub[i]);
       printf("s0 %f \n", cosangleparms_h->s0[i]);
    }
     */
    /*	printf("torsion\n");
       for (int i =0; i< torsionListSize; i++){
          printf("i %i\n", i);
          printf("kchi %f \n", torsionparms_h->kchi[i]);
          printf("delta %f \n", torsionparms_h->delta[i]);
       }*/
    /*
         printf("bpair\n");
    for (int i =0; i< bpairListSize; i++){
       printf("i %i\n", i);
       printf("sigma %f \n", bpairparms_h->sigma[i]);
       printf("eps %f \n", bpairparms_h->eps[i]);
       printf("shift %f \n", bpairparms_h->shift[i]);
    }
     */
    //assign residue type ids to each particle
    for (int i = 0; i < sys->nion; i++)
    {
        char *name = state->species[i]->name;
        RESI_CONN * resi_conn = findResiConnNew(charmmParms, name);
        //printf("particle %i %s %i \n", i, name, resi_conn->resID-1);
        bondparms_h->resID[i] = resi_conn->resID - 1;
    }

    printf("bl %i angle %i cosangle %i tors %i impr %i bpair %i\n", bondListSize, angleListSize, cosangleListSize, torsionListSize, improperListSize, bpairListSize);
}

void migrateBond(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    int bondListSize = parms->bondparms_h->bondListSize;

    gpu_memcpy_host2device(bond_g->numAtomsInResidue, bond_h->numAtomsInResidue, charmmParms->resiConnSize);
    gpu_memcpy_host2device(bond_g->numAtomsInResiduePrefixSum, bond_h->numAtomsInResiduePrefixSum, charmmParms->resiConnSize);
    gpu_memcpy_host2device(bond_g->resID, bond_h->resID, sys->nion);

    //bond related buffers
    gpu_memcpy_host2device((bond_g->bondListI), bond_h->bondListI, bondListSize);
    gpu_memcpy_host2device((bond_g->bondListJ), bond_h->bondListJ, bondListSize);
    gpu_memcpy_host2device((bond_g->resiBondListStarts), bond_h->resiBondListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((bond_g->resiBondListSizes), bond_h->resiBondListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((bond_g->atmBondListStarts), bond_h->atmBondListStarts, atomListSize);
    gpu_memcpy_host2device((bond_g->atmBondListSizes), bond_h->atmBondListSizes, atomListSize);
    gpu_memcpy_host2device((bond_g->kb), bond_h->kb, bondListSize);
    gpu_memcpy_host2device((bond_g->b0), bond_h->b0, bondListSize);


    gpu_memcpy_host2device(parms->gpu_bondparms, bond_g, 1);
}

void migrateAngle(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    //CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    //CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h;   //host pointer to host memory
    CharmmANGLEGPU_PARMS *angle_g = parms->gpu_angleparms_h; //host pointer to gpu memory 
    CharmmANGLEGPU_PARMS *angle_h = parms->angleparms_h; //host pointer to host memory


    int atomListSize = parms->bondparms_h->atomListSize;
    //int bondListSize = parms->bondparms_h->bondListSize;
    int angleListSize = parms->angleparms_h->angleListSize;

    //angle related buffers
    gpu_memcpy_host2device((angle_g->angleListI), angle_h->angleListI, angleListSize);
    gpu_memcpy_host2device((angle_g->angleListJ), angle_h->angleListJ, angleListSize);
    gpu_memcpy_host2device((angle_g->angleListK), angle_h->angleListK, angleListSize);
    gpu_memcpy_host2device((angle_g->resiAngleListStarts), angle_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((angle_g->resiAngleListSizes), angle_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((angle_g->atmAngleListStarts), angle_h->atmAngleListStarts, atomListSize);
    gpu_memcpy_host2device((angle_g->atmAngleListSizes), angle_h->atmAngleListSizes, atomListSize);
    gpu_memcpy_host2device((angle_g->ktheta), angle_h->ktheta, angleListSize);
    gpu_memcpy_host2device((angle_g->theta0), angle_h->theta0, angleListSize);
    gpu_memcpy_host2device((angle_g->kub), angle_h->kub, angleListSize);
    gpu_memcpy_host2device((angle_g->s0), angle_h->s0, angleListSize);

    gpu_memcpy_host2device(parms->gpu_angleparms, angle_g, 1);
}

void migrateCosAngle(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    //CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    //CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h;   //host pointer to host memory
    CharmmANGLEGPU_PARMS *cosangle_g = parms->gpu_cosangleparms_h; //host pointer to gpu memory 
    CharmmANGLEGPU_PARMS *cosangle_h = parms->cosangleparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    //int bondListSize = parms->bondparms_h->bondListSize;
    int cosangleListSize = parms->cosangleparms_h->angleListSize;

    //cos angle related buffers
    gpu_memcpy_host2device((cosangle_g->angleListI), cosangle_h->angleListI, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->angleListJ), cosangle_h->angleListJ, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->angleListK), cosangle_h->angleListK, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->resiAngleListStarts), cosangle_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((cosangle_g->resiAngleListSizes), cosangle_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((cosangle_g->atmAngleListStarts), cosangle_h->atmAngleListStarts, atomListSize);
    gpu_memcpy_host2device((cosangle_g->atmAngleListSizes), cosangle_h->atmAngleListSizes, atomListSize);
    gpu_memcpy_host2device((cosangle_g->ktheta), cosangle_h->ktheta, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->theta0), cosangle_h->theta0, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->kub), cosangle_h->kub, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->s0), cosangle_h->s0, cosangleListSize);

    gpu_memcpy_host2device(parms->gpu_cosangleparms, cosangle_g, 1);
}

void migrateRebAngle(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    CharmmANGLEGPU_PARMS *rebangle_g = parms->gpu_rebangleparms_h; //host pointer to gpu memory 
    CharmmANGLEGPU_PARMS *rebangle_h = parms->rebangleparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    //int bondListSize = parms->bondparms_h->bondListSize;
    int rebangleListSize = parms->rebangleparms_h->angleListSize;

    //Restricted bending angle related buffers
    gpu_memcpy_host2device((rebangle_g->angleListI), rebangle_h->angleListI, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->angleListJ), rebangle_h->angleListJ, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->angleListK), rebangle_h->angleListK, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->resiAngleListStarts), rebangle_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((rebangle_g->resiAngleListSizes), rebangle_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((rebangle_g->atmAngleListStarts), rebangle_h->atmAngleListStarts, atomListSize);
    gpu_memcpy_host2device((rebangle_g->atmAngleListSizes), rebangle_h->atmAngleListSizes, atomListSize);
    gpu_memcpy_host2device((rebangle_g->ktheta), rebangle_h->ktheta, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->theta0), rebangle_h->theta0, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->kub), rebangle_h->kub, rebangleListSize);
    gpu_memcpy_host2device((rebangle_g->s0), rebangle_h->s0, rebangleListSize);

    gpu_memcpy_host2device(parms->gpu_rebangleparms, rebangle_g, 1);
}

void migrateTorsion(SYSTEM * sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    //CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    //CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h;   //host pointer to host memory
    CharmmTORSIONGPU_PARMS *torsion_g = parms->gpu_torsionparms_h; //host pointer to gpu memory 
    CharmmTORSIONGPU_PARMS *torsion_h = parms->torsionparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    //int bondListSize = parms->bondparms_h->bondListSize;
    int torsionListSize = parms->torsionparms_h->torsionListSize;

    //torsion related buffers
    gpu_memcpy_host2device((torsion_g->torsionListI), torsion_h->torsionListI, torsionListSize);
    gpu_memcpy_host2device((torsion_g->torsionListJ), torsion_h->torsionListJ, torsionListSize);
    gpu_memcpy_host2device((torsion_g->torsionListK), torsion_h->torsionListK, torsionListSize);
    gpu_memcpy_host2device((torsion_g->torsionListL), torsion_h->torsionListL, torsionListSize);
    gpu_memcpy_host2device((torsion_g->resiTorsionListStarts), torsion_h->resiTorsionListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((torsion_g->resiTorsionListSizes), torsion_h->resiTorsionListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((torsion_g->atmTorsionListStarts), torsion_h->atmTorsionListStarts, atomListSize);
    gpu_memcpy_host2device((torsion_g->atmTorsionListSizes), torsion_h->atmTorsionListSizes, atomListSize);
    gpu_memcpy_host2device((torsion_g->kchi), torsion_h->kchi, torsionListSize);
    gpu_memcpy_host2device((torsion_g->delta), torsion_h->delta, torsionListSize);
    gpu_memcpy_host2device((torsion_g->n), torsion_h->n, torsionListSize);


    gpu_memcpy_host2device(parms->gpu_torsionparms, torsion_g, 1);

}

void migrateImproper(SYSTEM * sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    //CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    //CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h;   //host pointer to host memory
    CharmmIMPROPERGPU_PARMS *improper_g = parms->gpu_improperparms_h; //host pointer to gpu memory 
    CharmmIMPROPERGPU_PARMS *improper_h = parms->improperparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    //int bondListSize = parms->bondparms_h->bondListSize;
    int improperListSize = parms->improperparms_h->improperListSize;

    //improper related buffers
    gpu_memcpy_host2device((improper_g->improperListI), improper_h->improperListI, improperListSize);
    gpu_memcpy_host2device((improper_g->improperListJ), improper_h->improperListJ, improperListSize);
    gpu_memcpy_host2device((improper_g->improperListK), improper_h->improperListK, improperListSize);
    gpu_memcpy_host2device((improper_g->improperListL), improper_h->improperListL, improperListSize);
    gpu_memcpy_host2device((improper_g->resiImproperListStarts), improper_h->resiImproperListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((improper_g->resiImproperListSizes), improper_h->resiImproperListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((improper_g->atmImproperListStarts), improper_h->atmImproperListStarts, atomListSize);
    gpu_memcpy_host2device((improper_g->atmImproperListSizes), improper_h->atmImproperListSizes, atomListSize);
    gpu_memcpy_host2device((improper_g->kpsi), improper_h->kpsi, improperListSize);
    gpu_memcpy_host2device((improper_g->psi0), improper_h->psi0, improperListSize);


    gpu_memcpy_host2device(parms->gpu_improperparms, improper_g, 1);

}

void migrateBPair(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    //CHARMM_PARMS* charmmParms = parms->charmmParms;

    //CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    //CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h;   //host pointer to host memory
    CharmmBPAIRGPU_PARMS *bpair_g = parms->gpu_bpairparms_h; //host pointer to gpu memory 
    CharmmBPAIRGPU_PARMS *bpair_h = parms->bpairparms_h; //host pointer to host memory

    int atomListSize = parms->bondparms_h->atomListSize;
    int bpairListSize = parms->bpairparms_h->bpairListSize;

    //bpair related buffers
    gpu_memcpy_host2device((bpair_g->bpairListI), bpair_h->bpairListI, bpairListSize);
    gpu_memcpy_host2device((bpair_g->bpairListJ), bpair_h->bpairListJ, bpairListSize);
    gpu_memcpy_host2device((bpair_g->atmBpairListStarts), bpair_h->atmBpairListStarts, atomListSize);
    gpu_memcpy_host2device((bpair_g->atmBpairListSizes), bpair_h->atmBpairListSizes, atomListSize);
    gpu_memcpy_host2device((bpair_g->shift), bpair_h->shift, bpairListSize);
    gpu_memcpy_host2device((bpair_g->sigma), bpair_h->sigma, bpairListSize);
    gpu_memcpy_host2device((bpair_g->eps), bpair_h->eps, bpairListSize);

    gpu_memcpy_host2device(parms->gpu_bpairparms, bpair_g, 1);

}

void migrateResiCon(SYSTEM * sys, CHARMMPOT_PARMS *parms)
{
    migrateBond(sys, parms);
    migrateAngle(sys, parms);
    migrateCosAngle(sys, parms);
    migrateRebAngle(sys, parms);
    migrateTorsion(sys, parms);
    migrateImproper(sys, parms);
    migrateBPair(sys, parms);
}


// This function is deprecated
void migrateResiConOld(SYSTEM *sys, CHARMMPOT_PARMS *parms)
{
    //STATE * state = sys->collection->state;
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    CharmmBONDGPU_PARMS *bond_g = parms->gpu_bondparms_h; //host pointer to gpu memory
    CharmmBONDGPU_PARMS *bond_h = parms->bondparms_h; //host pointer to host memory
    CharmmANGLEGPU_PARMS *angle_g = parms->gpu_angleparms_h; //host pointer to gpu memory 
    CharmmANGLEGPU_PARMS *angle_h = parms->angleparms_h; //host pointer to host memory
    CharmmANGLEGPU_PARMS *cosangle_g = parms->gpu_cosangleparms_h; //host pointer to gpu memory 
    CharmmANGLEGPU_PARMS *cosangle_h = parms->cosangleparms_h; //host pointer to host memory


    int atomListSize = parms->bondparms_h->atomListSize;
    int bondListSize = parms->bondparms_h->bondListSize;
    int angleListSize = parms->angleparms_h->angleListSize;
    int cosangleListSize = parms->cosangleparms_h->angleListSize;

    gpu_memcpy_host2device(bond_g->numAtomsInResidue, bond_h->numAtomsInResidue, charmmParms->resiConnSize);
    gpu_memcpy_host2device(bond_g->numAtomsInResiduePrefixSum, bond_h->numAtomsInResiduePrefixSum, charmmParms->resiConnSize);
    gpu_memcpy_host2device(bond_g->resID, bond_h->resID, sys->nion);

    //bond related buffers
    gpu_memcpy_host2device((bond_g->bondListI), bond_h->bondListI, bondListSize);
    gpu_memcpy_host2device((bond_g->bondListJ), bond_h->bondListJ, bondListSize);
    gpu_memcpy_host2device((bond_g->resiBondListStarts), bond_h->resiBondListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((bond_g->resiBondListSizes), bond_h->resiBondListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((bond_g->atmBondListStarts), bond_h->atmBondListStarts, atomListSize);
    gpu_memcpy_host2device((bond_g->atmBondListSizes), bond_h->atmBondListSizes, atomListSize);
    gpu_memcpy_host2device((bond_g->kb), bond_h->kb, bondListSize);
    gpu_memcpy_host2device((bond_g->b0), bond_h->b0, bondListSize);


    //angle related buffers
    gpu_memcpy_host2device((angle_g->angleListI), angle_h->angleListI, angleListSize);
    gpu_memcpy_host2device((angle_g->angleListJ), angle_h->angleListJ, angleListSize);
    gpu_memcpy_host2device((angle_g->angleListK), angle_h->angleListK, angleListSize);
    gpu_memcpy_host2device((angle_g->resiAngleListStarts), angle_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((angle_g->resiAngleListSizes), angle_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((angle_g->atmAngleListStarts), angle_h->atmAngleListStarts, atomListSize);
    gpu_memcpy_host2device((angle_g->atmAngleListSizes), angle_h->atmAngleListSizes, atomListSize);
    gpu_memcpy_host2device((angle_g->ktheta), angle_h->ktheta, angleListSize);
    gpu_memcpy_host2device((angle_g->theta0), angle_h->theta0, angleListSize);
    gpu_memcpy_host2device((angle_g->kub), angle_h->kub, angleListSize);
    gpu_memcpy_host2device((angle_g->s0), angle_h->s0, angleListSize);


    //cos angle related buffers
    gpu_memcpy_host2device((cosangle_g->angleListI), cosangle_h->angleListI, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->angleListJ), cosangle_h->angleListJ, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->angleListK), cosangle_h->angleListK, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->resiAngleListStarts), cosangle_h->resiAngleListStarts, charmmParms->resiConnSize);
    gpu_memcpy_host2device((cosangle_g->resiAngleListSizes), cosangle_h->resiAngleListSizes, charmmParms->resiConnSize);
    gpu_memcpy_host2device((cosangle_g->atmAngleListStarts), cosangle_h->atmAngleListStarts, atomListSize);
    gpu_memcpy_host2device((cosangle_g->atmAngleListSizes), cosangle_h->atmAngleListSizes, atomListSize);
    gpu_memcpy_host2device((cosangle_g->ktheta), cosangle_h->ktheta, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->theta0), cosangle_h->theta0, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->kub), cosangle_h->kub, cosangleListSize);
    gpu_memcpy_host2device((cosangle_g->s0), cosangle_h->s0, cosangleListSize);


    //migrate helper structs
    gpu_memcpy_host2device(parms->gpu_bondparms, bond_g, 1);
    gpu_memcpy_host2device(parms->gpu_angleparms, angle_g, 1);
    gpu_memcpy_host2device(parms->gpu_cosangleparms, cosangle_g, 1);
}

void radix_sort_cub(int *r_backBondinp, int *r_backBond, gid_type* gids, gid_type* gidsSort, int nLocal, int nIon)
{
    //unsigned long long * gidsL = (unsigned long long*) gids;
    unsigned long long *buffer = NULL;
    unsigned long long * gidsSortL = (unsigned long long *) gidsSort;
    unsigned long long * gidsL = (unsigned long long *) gids;
    size_t temp_storage_bytes = 0;
    printf("starting sort\n");
    /*
    for (int i =0 ; i<nIon; i++)
    {
       printf("%i gid %lu \n", i, gids[i]);
    }
     */
    printf("end\n");
    CUDA_SAFE_CALL(cub::DeviceRadixSort::SortPairs(buffer, temp_storage_bytes, gidsL, gidsSortL, r_backBondinp, r_backBond, nIon);)
    CUDA_SAFE_CALL(cudaMalloc(&buffer, temp_storage_bytes);)
    CUDA_SAFE_CALL(cub::DeviceRadixSort::SortPairs(buffer, temp_storage_bytes, gidsL, gidsSortL, r_backBondinp, r_backBond, nIon);)
    cudaStreamSynchronize(0);

    //unsigned long long  * gidsC = (unsigned long long *) malloc(nIon*sizeof(unsigned long long));
    //cudaMemcpy(gidsC, gidsSortL, sizeof(unsigned long long )*nIon, cudaMemcpyDeviceToHost);
    /*
    for (int i =0 ; i<nIon; i++)
    {
       printf("%i gid %lu \n", i, gidsSortL[i]);
    }
     */
}



#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8

__device__ void preduceGPUO3b(THREE_MATRIX *h1, double *x, double *y, double *z)
{
    double *h = (double *) h1;
    double hxx, hyy, hzz;
    hxx = h[XX];
    hyy = h[YY];
    hzz = h[ZZ];

    if (*x > 0.5 * hxx)
    {
        *x += -hxx;
    }
    if (*x < -0.5 * hxx)
    {
        *x += hxx;
    }

    if (*y > 0.5 * hyy)
    {
        *y += -hyy;
    }
    if (*y < -0.5 * hyy)
    {
        *y += hyy;
    }

    if (*z > 0.5 * hzz)
    {
        *z += -hzz;
    }
    if (*z < -0.5 * hzz)
    {
        *z += hzz;
    }


}

//TODO
//we should eventually be able to use the same source for this and the cpu code

__device__ void bioVec(STATE* state, THREE_MATRIX *h1, int atom1, int atom2, THREE_VECTOR* vec)
{
    double bx = state->rx[atom1] - state->rx[atom2];
    double by = state->ry[atom1] - state->ry[atom2];
    double bz = state->rz[atom1] - state->rz[atom2];

    PBC(preduceGPUO3b(h1, &bx, &by, &bz);)

    vec->x = bx;
    vec->y = by;
    vec->z = bz;
}

__global__ void permuteParticlesWithGid(STATE * state, double *rxs, double *rys, double *rzs, double *qs, gid_type * gid, gid_type * gidS, int *r_backBond, int *resID, int *resIDs, int nIon)
{

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    // if (pid==1) printf("rx %p ry %p rz %p rxs %p rys %p rzs %p gids %p gid %p qs %p stateq %p \n", state->rx, state->ry, state->rz, rxs, rys, rzs, gidS, gid, qs, state->q); 
    int newID = r_backBond[pid];
    rxs[pid] = state->rx[newID];
    rys[pid] = state->ry[newID];
    rzs[pid] = state->rz[newID];
    gidS[pid] = gid[newID];
    qs[pid] = state->q[newID];
    //resIDs[pid] = resID[newID];

    return;
}


//__global__ void bondKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resid, int * r_back, double *e, GPUVIRIALS *sion, int * bondListI, int *bondListJ, int * bondListStarts, int * bondListSizes, int *atmBondListStarts, int * atmBondListSizes,int *numAtomsInResiduePrefixSum, double *b0, double *kb, THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void bondKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resid, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * bondListI, int *bondListJ, int * bondListStarts, int * bondListSizes, int *atmBondListStarts, int * atmBondListSizes, int *numAtomsInResiduePrefixSum, double *b0, double *kb, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int ii = r_back[pid];

    long long unsigned label = labelS[pid];
    int resID = resid[ii];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;

    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int bond_start_idx = atmBondListStarts[atm_idx];
    int bond_end_idx = atmBondListSizes[atm_idx] + bond_start_idx;

    //printf("ii %i resid %i start %i\n", ii,resID, bond_start_idx);
    if (bond_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((ii >= nLocal))
    {
        return;
    }
    //printf("ii %i bond_start %i\n", ii, bond_start_idx);
    double ix = rx[pid];
    double iy = ry[pid];
    double iz = rz[pid];
    double eb = 0;

    double ffx = 0;
    double ffy = 0;
    double ffz = 0;
    double fpxx = 0;
    double fpyy = 0;
    double fpzz = 0;
    double fpxy = 0;
    double fpxz = 0;
    double fpyz = 0;

    for (int bi = bond_start_idx; bi < bond_end_idx; bi++)
    {
        int indexJ = pid - bondListI[bi] + bondListJ[bi];
        double jx = rx[indexJ];
        double jy = ry[indexJ];
        double jz = rz[indexJ];

        double bx = ix - jx;
        double by = iy - jy;
        double bz = iz - jz;
        //double bx1=bx;
        //double by1=by;
        //double bz1=bz;
        PBC(preduceGPUO3b(h1, &bx, &by, &bz);)

            double kbb = kb[bi];
        double r2 = bx * bx + by * by + bz*bz;
        double rinv = rsqrt(r2);
        double r = r2*rinv;
        double bDelta = r - b0[bi];
        double kbbDelta = kbb*bDelta;
        double ee = kbbDelta*bDelta;
        double kforce = -2 * kbbDelta;

        //if(r<3.1) printf("pid=%d icoor= %f %f %f jcoor= %f %f %f b0=%f  r=%f e=%f\n", pid, ix, iy, iz, jx, jy, jz, b0[bi], r, ee);
        //if(pid<20) printf("pid=%d icoor= %f %f %f jcoor= %f %f %f b0=%f  r=%f \n", pid, ix, iy, iz, jx, jy, jz, b0[bi], r);

        double unit_x = bx*rinv;
        double unit_y = by*rinv;
        double unit_z = bz*rinv;
        double fxDelta = kforce*unit_x;
        double fyDelta = kforce*unit_y;
        double fzDelta = kforce*unit_z;
        ffx += fxDelta;
        ffy += fyDelta;
        ffz += fzDelta;
        eb += ee;
        int jj = r_back[indexJ];

        fpxx -= fxDelta*bx;
        fpyy -= fyDelta*by;
        fpzz -= fzDelta*bz;
        fpxy -= fxDelta*by;
        fpxz -= fxDelta*bz;
        fpyz -= fyDelta*bz;
        //atomicAdd(sion->xy+jj, -fpxy);
        //atomicAdd(sion->xz+jj, -fpxz);
        //atomicAdd(sion->yz+jj, -fpyz);


        atomicAdd(fx + jj, -fxDelta);
        atomicAdd(fy + jj, -fyDelta);
        atomicAdd(fz + jj, -fzDelta);
        //printf("ii %i j %i rx1 %f rx2 %f ee %f fx %f fy %f fz %f\n", ii, r_back[indexJ],ix,rx[indexJ],  ee, fxDelta, fyDelta, fzDelta);
    }
    e[ii] += eb;
    //atomicAdd(e+ii, eb);
    //if(eb>1.0 || eb<0) printf("pid = %d   ii= %d   eb= %f\n", pid, ii, eb);

    //sion->xx[ii]+= fpxx;
    //sion->yy[ii]+= fpyy;
    //sion->zz[ii]+= fpzz;
    //sion->xy[ii]+= fpxy;
    //sion->xz[ii]+= fpxz;
    //sion->yz[ii]+= fpyz;
    atomicAdd(sxx + ii, fpxx);
    atomicAdd(syy + ii, fpyy);
    atomicAdd(szz + ii, fpzz);
    /*
    sxx[ii]+= fpxx;
    syy[ii]+= fpyy;
    szz[ii]+= fpzz;
    sxy[ii]+= fpxy;
    sxz[ii]+= fpxz;
    syz[ii]+= fpyz;
     */
    atomicAdd(fx + ii, ffx);
    atomicAdd(fy + ii, ffy);
    atomicAdd(fz + ii, ffz);




    /*
      fx[ii]+=ffx;  
      fy[ii]+=ffy;  
      fz[ii]+=ffz;  
     */
    return;
}


//__global__ void angleKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e,GPUVIRIALS *sion, int * angleListI, int *angleListJ, int *angleListK, int * angleListStarts, int * angleListSizes, int *atmAngleListStarts, int * atmAngleListSizes,int *numAtomsInResiduePrefixSum, double *ktheta, double *theta0, double *kub, double *s0,THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void angleKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * angleListI, int *angleListJ, int *angleListK, int * angleListStarts, int * angleListSizes, int *atmAngleListStarts, int * atmAngleListSizes, int *numAtomsInResiduePrefixSum, double *ktheta, double *theta0, double *kub, double *s0, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int jj = r_back[pid];
    if (jj >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[jj];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;

    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int angle_start_idx = atmAngleListStarts[atm_idx];
    int angle_end_idx = atmAngleListSizes[atm_idx] + angle_start_idx;

    if (angle_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((jj >= nLocal))
    {
        return;
    }

    //printf("ii %i resid %i start %i end %i\n", jj,resID, angle_start_idx, angle_end_idx);

    //atom J is used for determining which atom in angle owns this force term
    double jx = rx[pid];
    double jy = ry[pid];
    double jz = rz[pid];
    double fxj = 0;
    double fyj = 0;
    double fzj = 0;

    double fpxx = 0;
    double fpyy = 0;
    double fpzz = 0; // double fpxy=0; double fpxz=0; double fpyz=0;
    double ea = 0;
    for (int bi = angle_start_idx; bi < angle_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        int indexI = pid - angleListJ[bi] + angleListI[bi];
        int indexK = pid - angleListJ[bi] + angleListK[bi];
        int ii = r_back[indexI];
        int kk = r_back[indexK];

        double ix = rx[indexI];
        double iy = ry[indexI];
        double iz = rz[indexI];

        double ijx = ix - jx;
        double ijy = iy - jy;
        double ijz = iz - jz;
        PBC(preduceGPUO3b(h1, &ijx, &ijy, &ijz);)
            double vec_ijx = ijx;
        double vec_ijy = ijy;
        double vec_ijz = ijz;

        double kx = rx[indexK];
        double ky = ry[indexK];
        double kz = rz[indexK];

        double kjx = kx - jx;
        double kjy = ky - jy;
        double kjz = kz - jz;

        PBC(preduceGPUO3b(h1, &kjx, &kjy, &kjz);)
            double vec_kjx = kjx;
        double vec_kjy = kjy;
        double vec_kjz = kjz;

        //calculate distance between particles I,J and J,K
        double ijr2 = ijx * ijx + ijy * ijy + ijz*ijz;
        double rinv_ij = 1 / sqrt(ijr2); // rsqrt(ijr2);  
        double ijr = sqrt(ijr2); //ijr2*rinv_ij;

        double kjr2 = kjx * kjx + kjy * kjy + kjz*kjz;
        double rinv_kj = 1 / sqrt(kjr2); //rsqrt(kjr2); 
        double kjr = sqrt(kjr2); // kjr2*rinv_kj;

        //normalize IJ and JK vectors to get unit vectors
        ijx = ijx*rinv_ij;
        ijy = ijy*rinv_ij;
        ijz = ijz*rinv_ij;

        kjx = kjx*rinv_kj;
        kjy = kjy*rinv_kj;
        kjz = kjz*rinv_kj;

        //hook's law
        double cosTheta = ijx * kjx + ijy * kjy + ijz*kjz; //dot product normalized(IJ) * normalized(JK) =cos(theta_ijk)
        double a = acos(cosTheta); //take the inverse cosine of the dot product to get theta of angle IJK=theta_ijk
        double aDelta = a - theta0[bi]; //difference between current angle and resting angle 
        double k_theta = ktheta[bi];
        double ee = k_theta * aDelta*aDelta;
        ea += ee; //square and multiply difference by angle "spring" constant


        double sinabs = sin(a);
        //if(sinabs<0) sinabs=-sinabs;
        double coef_i = 2 * k_theta * aDelta / (ijr * sinabs);
        double fxDeltaI = coef_i * (kjx - ijx * cosTheta);
        double fyDeltaI = coef_i * (kjy - ijy * cosTheta);
        double fzDeltaI = coef_i * (kjz - ijz * cosTheta);

        double coef_k = 2 * k_theta * aDelta / (kjr * sinabs);
        double fxDeltaK = coef_k * (ijx - kjx * cosTheta);
        double fyDeltaK = coef_k * (ijy - kjy * cosTheta);
        double fzDeltaK = coef_k * (ijz - kjz * cosTheta);

        //printf("ii %i jj %i ee %f fx %f fy %f fz %f rij %f costheta %f cos %f\n", ii, jj, ee, fxDeltaK, fyDeltaK, fzDeltaK,ijr, kx, cosTheta);
        atomicAdd(fx + ii, fxDeltaI);
        atomicAdd(fy + ii, fyDeltaI);
        atomicAdd(fz + ii, fzDeltaI);

        atomicAdd(fx + kk, fxDeltaK);
        atomicAdd(fy + kk, fyDeltaK);
        atomicAdd(fz + kk, fzDeltaK);

        fxj -= (fxDeltaI + fxDeltaK);
        fyj -= (fyDeltaI + fyDeltaK);
        fzj -= (fzDeltaI + fzDeltaK);


        // virials
        fpxx -= (fxDeltaI * vec_ijx + fxDeltaK * vec_kjx);
        fpyy -= (fyDeltaI * vec_ijy + fyDeltaK * vec_kjy);
        fpzz -= (fzDeltaI * vec_ijz + fzDeltaK * vec_kjz);
        // fpxy = (fxDeltaI*vec_ijy+fxDeltaK*vec_kjy);
        // fpxz = (fxDeltaI*vec_ijz+fxDeltaK*vec_kjz);
        // fpyz = (fyDeltaI*vec_ijz+fyDeltaK*vec_kjz);
        //  int j = pid;
        //aatomicAdd(sion->xz+j, -fpxz);
        //atomicAdd(sion->yz+j, -fpyz);
        //atomicAdd(sion->xy+j, -fpxy);


    }
    //sion->xx[jj]+=fpxx;
    //sion->yy[jj]+=fpyy;
    //sion->zz[jj]+=fpzz;
    //sion->xy[jj]+=fpxy;
    //sion->xz[jj]+=fpxz;
    //sion->yz[jj]+=fpyz;

    atomicAdd(sxx + jj, fpxx);
    atomicAdd(syy + jj, fpyy);
    atomicAdd(szz + jj, fpzz);
    /*
               sxx[jj]+=fpxx;
               syy[jj]+=fpyy;
               szz[jj]+=fpzz;
               sxy[jj]+=fpxy;
               sxz[jj]+=fpxz;
               syz[jj]+=fpyz;
     */
    atomicAdd(fx + jj, fxj);
    atomicAdd(fy + jj, fyj);
    atomicAdd(fz + jj, fzj);


    e[jj] += ea;

    return;
}



//__global__ void angleCosineKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, GPUVIRIALS *sion, int * angleListI, int *angleListJ, int *angleListK, int * angleListStarts, int * angleListSizes, int *atmAngleListStarts, int * atmAngleListSizes,int *numAtomsInResiduePrefixSum, double *ktheta, double *theta0, double *kub, double *s0, THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void angleCosineKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * angleListI, int *angleListJ, int *angleListK, int * angleListStarts, int * angleListSizes, int *atmAngleListStarts, int * atmAngleListSizes, int *numAtomsInResiduePrefixSum, double *ktheta, double *theta0, double *kub, double *s0, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int jj = r_back[pid];
    if (jj >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[jj];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;

    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int angle_start_idx = atmAngleListStarts[atm_idx];
    int angle_end_idx = atmAngleListSizes[atm_idx] + angle_start_idx;

    //printf("ii %i resid %i start %i end %i\n", jj,resID, angle_start_idx, angle_end_idx);
    if (angle_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((jj >= nLocal))
    {
        return;
    }


    //atom J is used for determining which atom in angle owns this force term
    double jx = rx[pid];
    double jy = ry[pid];
    double jz = rz[pid];
    double fxj = 0;
    double fyj = 0;
    double fzj = 0;

    double ea = 0;
    for (int bi = angle_start_idx; bi < angle_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        double fpxx = 0;
        double fpyy = 0;
        double fpzz = 0;
        double fpxy = 0;
        double fpxz = 0;
        double fpyz = 0;
        int indexI = pid - angleListJ[bi] + angleListI[bi];
        int indexK = pid - angleListJ[bi] + angleListK[bi];
        int ii = r_back[indexI];
        int kk = r_back[indexK];

        double ix = rx[indexI];
        double iy = ry[indexI];
        double iz = rz[indexI];

        double ijx = ix - jx;
        double ijy = iy - jy;
        double ijz = iz - jz;

        PBC(preduceGPUO3b(h1, &ijx, &ijy, &ijz);)
            double vec_ijx = ijx;
        double vec_ijy = ijy;
        double vec_ijz = ijz;

        double kx = rx[indexK];
        double ky = ry[indexK];
        double kz = rz[indexK];

        double kjx = kx - jx;
        double kjy = ky - jy;
        double kjz = kz - jz;

        PBC(preduceGPUO3b(h1, &kjx, &kjy, &kjz);)
            double vec_kjx = kjx;
        double vec_kjy = kjy;
        double vec_kjz = kjz;

        //calculate distance between particles I,J and J,K
        double ijr2 = ijx * ijx + ijy * ijy + ijz*ijz;
        double rinv_ij = 1 / sqrt(ijr2); // rsqrt(ijr2);  
        double ijr = sqrt(ijr2); //ijr2*rinv_ij;

        double kjr2 = kjx * kjx + kjy * kjy + kjz*kjz;
        double rinv_kj = 1 / sqrt(kjr2); //rsqrt(kjr2); 
        double kjr = sqrt(kjr2); // kjr2*rinv_kj;

        //normalize IJ and JK vectors to get unit vectors
        ijx = ijx*rinv_ij;
        ijy = ijy*rinv_ij;
        ijz = ijz*rinv_ij;

        kjx = kjx*rinv_kj;
        kjy = kjy*rinv_kj;
        kjz = kjz*rinv_kj;

        //hook's law
        double cosTheta = ijx * kjx + ijy * kjy + ijz*kjz; //dot product normalized(IJ) * normalized(JK) =cos(theta_ijk)
        double a = cosTheta;
        double aDelta = a - theta0[bi]; //difference between current cosine and resting cosine 
        double k_theta = ktheta[bi];
        double ee = k_theta * aDelta*aDelta;
        ea += ee; //square and multiply difference by "spring" constant


        double sinabs = sin(a);
        //if(sinabs<0) sinabs=-sinabs;
        double coef_i = -2 * k_theta * aDelta / (ijr);
        double fxDeltaI = coef_i * (kjx - ijx * cosTheta);
        double fyDeltaI = coef_i * (kjy - ijy * cosTheta);
        double fzDeltaI = coef_i * (kjz - ijz * cosTheta);

        double coef_k = -2 * k_theta * aDelta / (kjr);
        double fxDeltaK = coef_k * (ijx - kjx * cosTheta);
        double fyDeltaK = coef_k * (ijy - kjy * cosTheta);
        double fzDeltaK = coef_k * (ijz - kjz * cosTheta);

        //printf("ii %i jj %i ee %f fx %f fy %f fz %f rij %f costheta %f cos %f\n", ii, jj, ee, fxDeltaK, fyDeltaK, fzDeltaK,ijr, kx, cosTheta);
        atomicAdd(fx + ii, fxDeltaI);
        atomicAdd(fy + ii, fyDeltaI);
        atomicAdd(fz + ii, fzDeltaI);

        atomicAdd(fx + kk, fxDeltaK);
        atomicAdd(fy + kk, fyDeltaK);
        atomicAdd(fz + kk, fzDeltaK);

        fxj -= (fxDeltaI + fxDeltaK);
        fyj -= (fyDeltaI + fyDeltaK);
        fzj -= (fzDeltaI + fzDeltaK);

        // Pressure

        fpxx -= (fxDeltaI * vec_ijx + fxDeltaK * vec_kjx);
        fpxy -= (fxDeltaI * vec_ijy + fxDeltaK * vec_kjy);
        fpxz -= (fxDeltaI * vec_ijz + fxDeltaK * vec_kjz);
        fpyy -= (fyDeltaI * vec_ijy + fyDeltaK * vec_kjy);
        fpyz -= (fyDeltaI * vec_ijz + fyDeltaK * vec_kjz);
        fpzz -= (fzDeltaI * vec_ijz + fzDeltaK * vec_kjz);

        //atomicAdd(sion->xx+jj, fpxx);
        //atomicAdd(sion->yy+jj, fpyy);
        //atomicAdd(sion->zz+jj, fpzz);
        //atomicAdd(sion->xy+jj, fpxy);
        //atomicAdd(sion->xz+jj, fpxz);
        //atomicAdd(sion->yz+jj, fpyz);

        atomicAdd(sxx + jj, fpxx);
        atomicAdd(syy + jj, fpyy);
        atomicAdd(szz + jj, fpzz);
        atomicAdd(sxy + jj, fpxy);
        atomicAdd(sxz + jj, fpxz);
        atomicAdd(syz + jj, fpyz);


    }

    atomicAdd(fx + jj, fxj);
    atomicAdd(fy + jj, fyj);
    atomicAdd(fz + jj, fzj);
    e[jj] += ea;


    return;
}

__global__ void angleRebKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * angleListI, int *angleListJ, int *angleListK, int * angleListStarts, int * angleListSizes, int *atmAngleListStarts, int * atmAngleListSizes, int *numAtomsInResiduePrefixSum, double *ktheta, double *theta0, double *kub, double *s0, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int jj = r_back[pid];
    if (jj >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[jj];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;

    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int angle_start_idx = atmAngleListStarts[atm_idx];
    int angle_end_idx = atmAngleListSizes[atm_idx] + angle_start_idx;

    //printf("ii %i resid %i start %i end %i\n", jj,resID, angle_start_idx, angle_end_idx);
    if (angle_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((jj >= nLocal))
    {
        return;
    }


    //atom J is used for determining which atom in angle owns this force term
    double jx = rx[pid];
    double jy = ry[pid];
    double jz = rz[pid];
    double fxj = 0;
    double fyj = 0;
    double fzj = 0;

    double ea = 0;
    for (int bi = angle_start_idx; bi < angle_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        double fpxx = 0;
        double fpyy = 0;
        double fpzz = 0;
        double fpxy = 0;
        double fpxz = 0;
        double fpyz = 0;
        int indexI = pid - angleListJ[bi] + angleListI[bi];
        int indexK = pid - angleListJ[bi] + angleListK[bi];
        int ii = r_back[indexI];
        int kk = r_back[indexK];

        double ix = rx[indexI];
        double iy = ry[indexI];
        double iz = rz[indexI];

        double ijx = ix - jx;
        double ijy = iy - jy;
        double ijz = iz - jz;

        PBC(preduceGPUO3b(h1, &ijx, &ijy, &ijz);)
            double vec_ijx = ijx;
        double vec_ijy = ijy;
        double vec_ijz = ijz;

        double kx = rx[indexK];
        double ky = ry[indexK];
        double kz = rz[indexK];

        double kjx = kx - jx;
        double kjy = ky - jy;
        double kjz = kz - jz;

        PBC(preduceGPUO3b(h1, &kjx, &kjy, &kjz);)
            double vec_kjx = kjx;
        double vec_kjy = kjy;
        double vec_kjz = kjz;

        //calculate distance between particles I,J and J,K
        double ijr2 = ijx * ijx + ijy * ijy + ijz*ijz;
        double rinv_ij = 1 / sqrt(ijr2); // rsqrt(ijr2);  
        double ijr = sqrt(ijr2); //ijr2*rinv_ij;

        double kjr2 = kjx * kjx + kjy * kjy + kjz*kjz;
        double rinv_kj = 1 / sqrt(kjr2); //rsqrt(kjr2); 
        double kjr = sqrt(kjr2); // kjr2*rinv_kj;

        //normalize IJ and JK vectors to get unit vectors
        ijx = ijx*rinv_ij;
        ijy = ijy*rinv_ij;
        ijz = ijz*rinv_ij;

        kjx = kjx*rinv_kj;
        kjy = kjy*rinv_kj;
        kjz = kjz*rinv_kj;

        //hook's law
        double cosTheta = ijx * kjx + ijy * kjy + ijz*kjz; //dot product normalized(IJ) * normalized(JK) =cos(theta_ijk)
        double sinThetaSq = 1-cosTheta*cosTheta;
        
        double costTheta0=theta0[bi];
        double aDelta = cosTheta - costTheta0; //difference between current cosine and resting cosine 
        double k_theta = ktheta[bi];
        double ee = k_theta * aDelta*aDelta/sinThetaSq;
        ea += ee; //square and multiply difference by "spring" constant

        double ceof_reb= -2 * k_theta*aDelta*(1-cosTheta*costTheta0)/(sinThetaSq*sinThetaSq);
        
        double coef_i = ceof_reb / (ijr);
        double fxDeltaI = coef_i * (kjx - ijx * cosTheta);
        double fyDeltaI = coef_i * (kjy - ijy * cosTheta);
        double fzDeltaI = coef_i * (kjz - ijz * cosTheta);

        double coef_k = ceof_reb / (kjr);
        double fxDeltaK = coef_k * (ijx - kjx * cosTheta);
        double fyDeltaK = coef_k * (ijy - kjy * cosTheta);
        double fzDeltaK = coef_k * (ijz - kjz * cosTheta);

        //printf("ii %i jj %i ee %f fx %f fy %f fz %f rij %f costheta %f cos %f\n", ii, jj, ee, fxDeltaK, fyDeltaK, fzDeltaK,ijr, kx, cosTheta);
        atomicAdd(fx + ii, fxDeltaI);
        atomicAdd(fy + ii, fyDeltaI);
        atomicAdd(fz + ii, fzDeltaI);

        atomicAdd(fx + kk, fxDeltaK);
        atomicAdd(fy + kk, fyDeltaK);
        atomicAdd(fz + kk, fzDeltaK);

        fxj -= (fxDeltaI + fxDeltaK);
        fyj -= (fyDeltaI + fyDeltaK);
        fzj -= (fzDeltaI + fzDeltaK);

        // Pressure

        fpxx -= (fxDeltaI * vec_ijx + fxDeltaK * vec_kjx);
        fpxy -= (fxDeltaI * vec_ijy + fxDeltaK * vec_kjy);
        fpxz -= (fxDeltaI * vec_ijz + fxDeltaK * vec_kjz);
        fpyy -= (fyDeltaI * vec_ijy + fyDeltaK * vec_kjy);
        fpyz -= (fyDeltaI * vec_ijz + fyDeltaK * vec_kjz);
        fpzz -= (fzDeltaI * vec_ijz + fzDeltaK * vec_kjz);

        //atomicAdd(sion->xx+jj, fpxx);
        //atomicAdd(sion->yy+jj, fpyy);
        //atomicAdd(sion->zz+jj, fpzz);
        //atomicAdd(sion->xy+jj, fpxy);
        //atomicAdd(sion->xz+jj, fpxz);
        //atomicAdd(sion->yz+jj, fpyz);

        atomicAdd(sxx + jj, fpxx);
        atomicAdd(syy + jj, fpyy);
        atomicAdd(szz + jj, fpzz);
        atomicAdd(sxy + jj, fpxy);
        atomicAdd(sxz + jj, fpxz);
        atomicAdd(syz + jj, fpyz);


    }

    atomicAdd(fx + jj, fxj);
    atomicAdd(fy + jj, fyj);
    atomicAdd(fz + jj, fzj);
    e[jj] += ea;


    return;
}

__device__ double bioNorm(THREE_VECTOR vec)
{
    return (sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z));
}

__device__ THREE_VECTOR cross_g(THREE_VECTOR*a, THREE_VECTOR*b)
{
    THREE_VECTOR c;
    c.x = a->y * b->z - a->z * b->y;
    c.y = a->z * b->x - a->x * b->z;
    c.z = a->x * b->y - a->y * b->x;
    return c;
}

__device__ double dot1_g(THREE_VECTOR a, THREE_VECTOR b)
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

__device__ THREE_VECTOR bioVecAdd(double scale1, THREE_VECTOR vec1, double scale2, THREE_VECTOR vec2)
{
    THREE_VECTOR vec;
    vec.x = scale1 * vec1.x + scale2 * vec2.x;
    vec.y = scale1 * vec1.y + scale2 * vec2.y;
    vec.z = scale1 * vec1.z + scale2 * vec2.z;
    return vec;
}

__device__ THREE_VECTOR bioScale(THREE_VECTOR vec, double scale)
{
    vec.x *= scale;
    vec.y *= scale;
    vec.z *= scale;

    return vec;
}

__device__ int bioDihedralFast(STATE* state, THREE_MATRIX *h1, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                               THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL, THREE_SMATRIX *sion)
{
    double eps = 1e-12;
    THREE_VECTOR v_ij;
    bioVec(state, h1, indexI, indexJ, &v_ij);
    THREE_VECTOR v_jk;
    bioVec(state, h1, indexJ, indexK, &v_jk);
    THREE_VECTOR v_kl;
    bioVec(state, h1, indexK, indexL, &v_kl);

    double a2 = VSQ(v_ij);
    double b2 = VSQ(v_jk);
    double c2 = VSQ(v_kl);
    double ab = DOT(v_ij, v_jk);
    double bc = DOT(v_jk, v_kl);
    double ac = DOT(v_ij, v_kl);
    double f = ab * bc - ac*b2;
    double g1 = a2 * b2 - SQ(ab) + eps;
    double g2 = b2 * c2 - SQ(bc) + eps;
    double y = 1.0 / sqrt(g1 * g2);
    double x = y*f;
    double xab = y * bc + x / g1*ab;
    double xbc = y * ab + x / g2*bc;
    double xac = -y*b2;
    double xaa = -0.5 * x * b2 / g1;
    double xcc = -0.5 * x * b2 / g2;
    double xbb = -y * ac - 0.5 * x * (a2 / g1 + c2 / g2);
    THREE_VECTOR ca, cb, cc;
    VSVOP(ca, =, xab, *, v_jk);
    VSVOP(ca, +=, xac, *, v_kl);
    VSVOP(ca, +=, (2 * xaa), *, v_ij);

    VSVOP(cb, =, xab, *, v_ij);
    VSVOP(cb, +=, xbc, *, v_kl);
    VSVOP(cb, +=, (2 * xbb), *, v_jk);

    VSVOP(cc, =, xbc, *, v_jk);
    VSVOP(cc, +=, xac, *, v_ij);
    VSVOP(cc, +=, (2 * xcc), *, v_kl);

    //Decide the sign of X angle
    THREE_VECTOR m;
    CROSS(m, v_ij, v_jk);
    THREE_VECTOR n;
    CROSS(n, v_jk, v_kl);
    THREE_VECTOR mXn;
    CROSS(mXn, m, n); // Flip the conventional in Bekker 1995 JCC.
    double signnum = DOT(v_jk, mXn);
    double sign = (signnum < 0.0) ? -1.0 : 1.0;
    x = max(min(x, 1.0), -1.0);
    *angleX = sign * acos(x);
    *sinXptr = sin(*angleX);

    (*dcosXdrI) = ca;
    VOP2((*dcosXdrJ), =, cb, -, ca);
    VOP2((*dcosXdrK), =, cc, -, cb);
    VSVOP((*dcosXdrL), =, -1, *, cc);

    sion->xx = -1 * (ca.x * v_ij.x + cb.x * v_jk.x + cc.x * v_kl.x);
    sion->xy = -1 * (ca.x * v_ij.y + cb.x * v_jk.y + cc.x * v_kl.y);
    sion->xz = -1 * (ca.x * v_ij.z + cb.x * v_jk.z + cc.x * v_kl.z);
    sion->yy = -1 * (ca.y * v_ij.y + cb.y * v_jk.y + cc.y * v_kl.y);
    sion->yz = -1 * (ca.y * v_ij.z + cb.y * v_jk.z + cc.y * v_kl.z);
    sion->zz = -1 * (ca.z * v_ij.z + cb.z * v_jk.z + cc.z * v_kl.z);

    return 0;
}

__device__ int bioDihedral(STATE* state, THREE_MATRIX *h1, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                           THREE_VECTOR* v_ij, THREE_VECTOR* v_kj, THREE_VECTOR* v_kl,
                           THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL)
{

    //THREE_VECTOR v_ij;
    bioVec(state, h1, indexI, indexJ, v_ij);
    //THREE_VECTOR v_kj;
    bioVec(state, h1, indexK, indexJ, v_kj);
    //THREE_VECTOR v_kl;
    bioVec(state, h1, indexK, indexL, v_kl);

    THREE_VECTOR v_A = cross_g(v_ij, v_kj);
    THREE_VECTOR v_B = cross_g(v_kj, v_kl);
    THREE_VECTOR v_C = cross_g(&v_A, v_kj);

    double v1_A = bioNorm(v_A);
    double v1_B = bioNorm(v_B);
    double v1_C = bioNorm(v_C);
    //printf("norms i %i j %i k %i l %i v1_A %f v1_b %f v1_C %f\n", indexI, indexJ, indexK, indexL, v1_A, v1_B, v1_C); 
    if (v1_A < 1.0e-7)
    {
        printf("Warning: A vector rij X rkj is too small\n");
    }
    if (v1_B < 1.0e-7)
    {
        printf("Warning: B vector rkj X rkl is too small\n");
    }
    if (v1_C < 1.0e-7)
    {
        printf("Warning: C vector rkj X A is too small\n");
    }


    double dotCB = dot1_g(v_C, v_B);
    double sinX = dotCB / (v1_C * v1_B);
    *sinXptr = sinX;

    double dotAB = dot1_g(v_A, v_B);
    double cosX = dotAB / (v1_A * v1_B);

    /*
        THREE_VECTOR cross_gAB=cross_g(&v_A, &v_B);
        double dotkjAB=dot1(*v_kj, cross_gAB);
    
        double sign = 1.0;
        if(dotkjAB<0) sign = -1.0;
    
     *angleX=sign*acos(cosX);
     */
    if (cosX > 0)
    {
        *angleX = atan(sinX / cosX);
    }
    else if (sinX >= 0 && cosX < 0)
    {
        *angleX = atan(sinX / cosX) + M_PI;
    }
    else if (sinX < 0 && cosX < 0)
    {
        *angleX = atan(sinX / cosX) - M_PI;
    }
    else if (sinX > 0 && cosX == 0)
    {
        *angleX = M_PI / 2;
    }
    else if (sinX < 0 && cosX == 0)
    {
        *angleX = -M_PI / 2;
    }
    else if (sinX == 0 && cosX == 0)
    {
        //normally won't reach here.
        printf("ERROR: undefined angle\n");
        *angleX = 0;
    }


    double diffsinX = fabs(sinX - sin(*angleX));
    if (diffsinX > 1e-4)
    {
        printf("ERROR: X angle %f is not consist with sinX, %f \n", *angleX, sinX);
    }

    double coef = 1 / (v1_A * v1_B);
    double coef_a = -cosX / (v1_A * v1_A);
    double coef_b = -cosX / (v1_B * v1_B);
    //      1   B     A               1  
    // a = ---(--- - --- cosX) * (- ----)
    //     |A| |B|   |A|            sinX
    // b are similar. (-1/sinX is from dX/dcosX)
    THREE_VECTOR v_a = bioVecAdd(coef, v_B, coef_a, v_A);
    THREE_VECTOR v_b = bioVecAdd(coef, v_A, coef_b, v_B);

    // dcosX/drI = - r_kj X a
    (*dcosXdrI) = cross_g(&v_a, v_kj);

    // dcosX/drL = - r_kj X b   
    (*dcosXdrL) = cross_g(&v_b, v_kj);

    // dcosX/drJ = - r_ik X a + r_kl X b
    THREE_VECTOR v_ik;
    bioVec(state, h1, indexI, indexK, &v_ik);
    THREE_VECTOR va_ik = cross_g(&v_ik, &v_a);
    THREE_VECTOR vb_kl = cross_g(v_kl, &v_b);
    (*dcosXdrJ) = bioVecAdd(-1, va_ik, 1, vb_kl);

    // dcosX/drK = - r_jl X b + r_ij X a
    THREE_VECTOR v_jl;
    bioVec(state, h1, indexJ, indexL, &v_jl);
    THREE_VECTOR vb_jl = cross_g(&v_jl, &v_b);
    THREE_VECTOR va_ij = cross_g(v_ij, &v_a);
    (*dcosXdrK) = bioVecAdd(-1, vb_jl, 1, va_ij);

    return 0;


}

__global__ void improperKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, GPUVIRIALS *sion, int * improperListI, int *improperListJ, int *improperListK, int *improperListL, int * improperListStarts, int * improperListSizes, int *atmImproperListStarts, int * atmImproperListSizes, int *numAtomsInResiduePrefixSum, double *kpsi_, double *psi0_, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int ii = r_back[pid];
    if (ii >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[ii];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;


    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int improper_start_idx = atmImproperListStarts[atm_idx];
    int improper_end_idx = atmImproperListSizes[atm_idx] + improper_start_idx;

    if (improper_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((ii >= nLocal))
    {
        return;
    }

    STATE gpustate;
    STATE * gstate = &gpustate;
    gstate->rx = rx;
    gstate->ry = ry;
    gstate->rz = rz;

    //atom J is used for determining which atom in angle owns this force term
    int indexI = pid;

    double fpxx = 0;
    double fpyy = 0;
    double fpzz = 0;
    double fpxy = 0;
    double fpxz = 0;
    double fpyz = 0;
    double eimprtot = 0;
    double PI2=2*M_PI;
    double PI_1=-1*M_PI;
    //double eimpr =0;
    for (int bi = improper_start_idx; bi < improper_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        //int indexI = pid - atm_idx + improperListI[bi];
        int indexJ = pid - improperListI[bi] + improperListJ[bi];
        int indexK = pid - improperListI[bi] + improperListK[bi];
        int indexL = pid - improperListI[bi] + improperListL[bi];
        //int ii = r_back[indexI];
        int jj = r_back[indexJ];
        int kk = r_back[indexK];
        int ll = r_back[indexL];
        // if (pid<=20)
        //printf("pid %i  atmid %i i %i j %i k %i l %i ii %i jj %i kk %i ll %i \n\t I%i J%i K%i L %i size %i\n", pid,atm_idx, indexI, indexJ, indexK, indexL, ii, jj, kk, ll, improperListI[bi], improperListJ[bi], improperListK[bi], improperListL[bi], atmImproperListSizes[atm_idx]);   
        //Energy
        double impr;
        double sinX;
        THREE_VECTOR v_ij;
        THREE_VECTOR v_kj;
        THREE_VECTOR v_kl;

        THREE_VECTOR dTdrI;
        THREE_VECTOR dTdrJ;
        THREE_VECTOR dTdrK;
        THREE_VECTOR dTdrL;


        if (bioDihedral(gstate, h1, indexI, indexJ, indexK, indexL, &impr, &sinX,
                        &v_ij, &v_kj, &v_kl, &dTdrI, &dTdrJ, &dTdrK, &dTdrL) < 0)
        {
            continue; // skip the torsion with angle approach to 0.
            // also skip the force calculation because a & b approach to 0;
        }

        double psi0 = psi0_[bi];
        double kpsi = kpsi_[bi];
        double imprDelta = impr - psi0;
        if(imprDelta<PI_1){
            imprDelta=imprDelta+PI2;
        }else if (imprDelta> M_PI){
            imprDelta=imprDelta-PI2;
        }

        //double imprDelta=impr-imprPtr->psi0/RAD2DEG;
        double eimpr = kpsi * imprDelta*imprDelta;
        eimprtot = eimprtot + eimpr;
        if (eimpr > 0.1)
        {
            // printf("I=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            // printf("Improper =%f Improper angle=%f kpsi=%f psi0=%f Improper Energy=%f\n",
            //		  impr, impr*RAD2DEG, kpsi, psi0*RAD2DEG, eimpr);
        }
        //Force
        //Improper is the permutation of Torsion.

        double absX = sinX;
        if (sinX < 0) absX = -sinX;

        double kimpr;

        if (absX > FLOAT_EPS)
        {
            kimpr = -2 * kpsi * imprDelta / sinX;
        }
        else
        {
            if (psi0 > NEAR_ZERO_ANGLE)
            {
                //  printf("Warning: imprPtr->kpsi doesn't equal to 0");
            }
            double impr2 = impr*impr;
            double impr4 = impr2*impr2;
            double impr6 = impr4*impr2;
            double impr8 = impr4*impr4;
            double impr10 = impr8*impr2;
            kimpr = -2 * kpsi / (1 - impr2 / 6 + impr4 / 120 - impr6 / 5040 + impr8 / 362880 - impr10 / 39916800);
        }


        kimpr = -kimpr; // Flip the sign for bioVec

        THREE_VECTOR FdTdRi = bioScale(dTdrI, kimpr);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, kimpr);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, kimpr);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, kimpr);

        /*
        fpxx += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.x+(FdTdRk.x-FdTdRj.x)*v_kj.x+(FdTdRk.x-FdTdRl.x)*v_kl.x);
        fpxy += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.y+(FdTdRk.x-FdTdRj.x)*v_kj.y+(FdTdRk.x-FdTdRl.x)*v_kl.y);
        fpxz += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.z+(FdTdRk.x-FdTdRj.x)*v_kj.z+(FdTdRk.x-FdTdRl.x)*v_kl.z); 
        fpyy += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.y+(FdTdRk.y-FdTdRj.y)*v_kj.y+(FdTdRk.y-FdTdRl.y)*v_kl.y);
        fpyz += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.z+(FdTdRk.y-FdTdRj.y)*v_kj.z+(FdTdRk.y-FdTdRl.y)*v_kl.z);
        fpzz += 0.5*((FdTdRi.z-FdTdRj.z)*v_ij.z+(FdTdRk.z-FdTdRj.z)*v_kj.z+(FdTdRk.z-FdTdRl.z)*v_kl.z);
         */

        atomicAdd(fx + ii, -FdTdRi.x);
        atomicAdd(fy + ii, -FdTdRi.y);
        atomicAdd(fz + ii, -FdTdRi.z);

        atomicAdd(fx + jj, -FdTdRj.x);
        atomicAdd(fy + jj, -FdTdRj.y);
        atomicAdd(fz + jj, -FdTdRj.z);

        atomicAdd(fx + kk, -FdTdRk.x);
        atomicAdd(fy + kk, -FdTdRk.y);
        atomicAdd(fz + kk, -FdTdRk.z);

        atomicAdd(fx + ll, -FdTdRl.x);
        atomicAdd(fy + ll, -FdTdRl.y);
        atomicAdd(fz + ll, -FdTdRl.z);

    }

    atomicAdd(sion->xx + ii, -fpxx);
    atomicAdd(sion->yy + ii, -fpyy);
    atomicAdd(sion->zz + ii, -fpzz);
    atomicAdd(sion->xy + ii, -fpxy);
    atomicAdd(sion->xz + ii, -fpxz);
    atomicAdd(sion->yz + ii, -fpyz);

    atomicAdd(e + ii, eimprtot);
}


//__global__ void improperKernelFast(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, GPUVIRIALS *sion, int * improperListI, int *improperListJ, int *improperListK, int *improperListL, int * improperListStarts, int * improperListSizes, int *atmImproperListStarts, int * atmImproperListSizes,int *numAtomsInResiduePrefixSum, double *kpsi_, double *psi0_, THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void improperKernelFast(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * improperListI, int *improperListJ, int *improperListK, int *improperListL, int * improperListStarts, int * improperListSizes, int *atmImproperListStarts, int * atmImproperListSizes, int *numAtomsInResiduePrefixSum, double *kpsi_, double *psi0_, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int ii = r_back[pid];
    if (ii >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[ii];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;


    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int improper_start_idx = atmImproperListStarts[atm_idx];
    int improper_end_idx = atmImproperListSizes[atm_idx] + improper_start_idx;

    if (improper_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((ii >= nLocal))
    {
        return;
    }

    STATE gpustate;
    STATE * gstate = &gpustate;
    gstate->rx = rx;
    gstate->ry = ry;
    gstate->rz = rz;

    //atom J is used for determining which atom in angle owns this force term
    int indexI = pid;

    //double fpxx=0; double fpyy=0; double fpzz=0; double fpxy=0; double fpxz=0; double fpyz=0;
    double eimprtot = 0;
    //double eimpr =0;
    for (int bi = improper_start_idx; bi < improper_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        //int indexI = pid - atm_idx + improperListI[bi];
        int indexJ = pid - improperListI[bi] + improperListJ[bi];
        int indexK = pid - improperListI[bi] + improperListK[bi];
        int indexL = pid - improperListI[bi] + improperListL[bi];
        //int ii = r_back[indexI];
        int jj = r_back[indexJ];
        int kk = r_back[indexK];
        int ll = r_back[indexL];

        double impr;
        double sinX;
        THREE_VECTOR dTdrI;
        THREE_VECTOR dTdrJ;
        THREE_VECTOR dTdrK;
        THREE_VECTOR dTdrL;
        THREE_SMATRIX ssion;
        bioDihedralFast(state, h1, indexI, indexJ, indexK, indexL, &impr, &sinX, &dTdrI, &dTdrJ, &dTdrK, &dTdrL, &ssion);


        double psi0 = psi0_[bi];
        double kpsi = kpsi_[bi];
        double imprDelta = impr - psi0;

        //double imprDelta=impr-imprPtr->psi0/RAD2DEG;
        double eimpr = kpsi * imprDelta*imprDelta;
        eimprtot = eimprtot + eimpr;
        if (eimpr > 0.1)
        {
            // printf("I=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            // printf("Improper =%f Improper angle=%f kpsi=%f psi0=%f Improper Energy=%f\n",
            //		  impr, impr*RAD2DEG, kpsi, psi0*RAD2DEG, eimpr);
        }
        //Force
        //Improper is the permutation of Torsion.

        double absX = sinX;
        if (sinX < 0) absX = -sinX;

        double kimpr;

        if (absX > FLOAT_EPS)
        {
            kimpr = -2 * kpsi * imprDelta / sinX;
        }
        else
        {
            if (psi0 > NEAR_ZERO_ANGLE)
            {
                //  printf("Warning: imprPtr->kpsi doesn't equal to 0");
            }
            double impr2 = impr*impr;
            double impr4 = impr2*impr2;
            double impr6 = impr4*impr2;
            double impr8 = impr4*impr4;
            double impr10 = impr8*impr2;
            kimpr = -2 * kpsi / (1 - impr2 / 6 + impr4 / 120 - impr6 / 5040 + impr8 / 362880 - impr10 / 39916800);
        }

        //kimpr=-kimpr; // Flip the sign for bioVec

        THREE_VECTOR FdTdRi = bioScale(dTdrI, kimpr);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, kimpr);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, kimpr);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, kimpr);

        atomicAdd(fx + ii, -FdTdRi.x);
        atomicAdd(fy + ii, -FdTdRi.y);
        atomicAdd(fz + ii, -FdTdRi.z);

        atomicAdd(fx + jj, -FdTdRj.x);
        atomicAdd(fy + jj, -FdTdRj.y);
        atomicAdd(fz + jj, -FdTdRj.z);

        atomicAdd(fx + kk, -FdTdRk.x);
        atomicAdd(fy + kk, -FdTdRk.y);
        atomicAdd(fz + kk, -FdTdRk.z);

        atomicAdd(fx + ll, -FdTdRl.x);
        atomicAdd(fy + ll, -FdTdRl.y);
        atomicAdd(fz + ll, -FdTdRl.z);

        //atomicAdd(sion->xx+ii, -ssion.xx*kimpr);
        //atomicAdd(sion->yy+ii, -ssion.yy*kimpr);
        //atomicAdd(sion->zz+ii, -ssion.zz*kimpr);
        //atomicAdd(sion->xy+ii, -ssion.xy*kimpr);
        //atomicAdd(sion->xz+ii, -ssion.xz*kimpr);
        //atomicAdd(sion->yz+ii, -ssion.yz*kimpr);

        atomicAdd(sxx + ii, -ssion.xx * kimpr);
        atomicAdd(syy + ii, -ssion.yy * kimpr);
        atomicAdd(szz + ii, -ssion.zz * kimpr);
        atomicAdd(sxy + ii, -ssion.xy * kimpr);
        atomicAdd(sxz + ii, -ssion.xz * kimpr);
        atomicAdd(syz + ii, -ssion.yz * kimpr);
    }


    atomicAdd(e + ii, eimprtot);
}

__global__ void torsionKernel(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, GPUVIRIALS *sion, int * torsionListI, int *torsionListJ, int *torsionListK, int *torsionListL, int * torsionListStarts, int * torsionListSizes, int *atmTorsionListStarts, int * atmTorsionListSizes, int *numAtomsInResiduePrefixSum, double *kchi1, double *delta1, int *n1, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int jj = r_back[pid];
    if (jj >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[jj];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;


    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int torsion_start_idx = atmTorsionListStarts[atm_idx];
    int torsion_end_idx = atmTorsionListSizes[atm_idx] + torsion_start_idx;

    //printf("jj %i atmidx %i  resid %i start %i end %i\n", jj,atm_idx,resID, torsion_start_idx, torsion_end_idx);
    if (torsion_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((jj >= nLocal))
    {
        return;
    }

    STATE gpustate;
    STATE * gstate = &gpustate;
    gstate->rx = rx;
    gstate->ry = ry;
    gstate->rz = rz;

    //atom J is used for determining which atom in angle owns this force term
    int indexJ = pid;

    //double fpxx=0; double fpyy=0; double fpzz=0; double fpxy=0; double fpxz=0; double fpyz=0;
    double etorstot = 0;
    //double et =0;
    for (int bi = torsion_start_idx; bi < torsion_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        int indexI = pid - torsionListJ[bi] + torsionListI[bi];
        int indexK = pid - torsionListJ[bi] + torsionListK[bi];
        int indexL = pid - torsionListJ[bi] + torsionListL[bi];
        int ii = r_back[indexI];
        int kk = r_back[indexK];
        int ll = r_back[indexL];

        //Energy
        double tors;
        double sinX;
        THREE_VECTOR v_ij;
        THREE_VECTOR v_kj;
        THREE_VECTOR v_kl;

        THREE_VECTOR dTdrI;
        THREE_VECTOR dTdrJ;
        THREE_VECTOR dTdrK;
        THREE_VECTOR dTdrL;

        double kchi = kchi1[bi];
        double delta = delta1[bi];
        int n = n1[bi];

        //printf("ii %i j %i k %i l %i rx1 %f rx2 %f kchi %f delta %f n %i\n", ii, r_back[indexJ],r_back[indexK], r_back[indexL],rx[indexI],rx[indexJ], kchi, delta,n);
        if (bioDihedral(gstate, h1, indexI, indexJ, indexK, indexL, &tors, &sinX,
                        &v_ij, &v_kj, &v_kl, &dTdrI, &dTdrJ, &dTdrK, &dTdrL) < 0)
        {
            continue; // skip the torsion with angle approach to 0.
            // also skip the force calculation because a & b approach to 0;
        }




        double etors = kchi * (1 + cos(n * tors - delta));
        etorstot = etorstot + etors;

        //printf("ii %i j %i rx1 %f rx2 %f tors %f etors %f \n", ii, r_back[indexJ],rx[indexI],rx[indexJ],  tors, etors);
        if (etors > 0.1 || etors<-0.1)
        {
            printf("I=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            //printf("Torsion =%f Torsion angle=%f kchi=%f n=%d delta=%f Torsion Energy=%f\n", 
            //    tors*RAD2DEG, tors2*RAD2DEG, torsPtr->kchi, torsPtr->n, torsPtr->delta*RAD2DEG, etors);  
        }
        //Force
        double absX = fabs(sinX);

        double ktors;

        if (absX > FLOAT_EPS)
        {
            ktors = kchi * n * sin(n * tors - delta) / sinX;
            //negative sign canceled out
        }
        else
        {
            double nX = n*tors;
            double nX2 = nX*nX;
            double nX4 = nX2*nX2;
            double nX6 = nX4*nX2;
            double nX8 = nX4*nX4;
            double nX10 = nX8*nX2;
            double X2 = tors*tors;
            double X4 = X2*X2;
            double X6 = X4*X2;
            double X8 = X4*X4;
            double X10 = X8*X2;
            double ratio = n * (1 - nX2 / 6 + nX4 / 120 - nX6 / 5040 + nX8 / 362880 - nX10 / 39916800) / (1 - X2 / 6 + X4 / 120 - X6 / 5040 + X8 / 362880 - X10 / 39916800);
            if (delta < NEAR_ZERO_ANGLE)
            { // for delta=0
                ktors = kchi * n*ratio;
            }
            else if (delta > NEAR_180_ANGLE)
            { // for delta=180
                ktors = -kchi * n*ratio;
                // sin(x-pi)=-sin(pi-x)=-sin(x)
            }
            else
            {
                //printf("Warning: torsPtr->delta doesn't equal to 0 or 180");
                ktors = kchi * n*ratio;
            }
        }


        ktors = -ktors; // Flip the sign for bioVec
        /*
        fpxx += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.x+(FdTdRk.x-FdTdRj.x)*v_kj.x+(FdTdRk.x-FdTdRl.x)*v_kl.x);
        fpxy += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.y+(FdTdRk.x-FdTdRj.x)*v_kj.y+(FdTdRk.x-FdTdRl.x)*v_kl.y);
        fpxz += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.z+(FdTdRk.x-FdTdRj.x)*v_kj.z+(FdTdRk.x-FdTdRl.x)*v_kl.z); 
        fpyy += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.y+(FdTdRk.y-FdTdRj.y)*v_kj.y+(FdTdRk.y-FdTdRl.y)*v_kl.y);
        fpyz += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.z+(FdTdRk.y-FdTdRj.y)*v_kj.z+(FdTdRk.y-FdTdRl.y)*v_kl.z);
        fpzz += 0.5*((FdTdRi.z-FdTdRj.z)*v_ij.z+(FdTdRk.z-FdTdRj.z)*v_kj.z+(FdTdRk.z-FdTdRl.z)*v_kl.z);
         */

        THREE_VECTOR FdTdRi = bioScale(dTdrI, ktors);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, ktors);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, ktors);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, ktors);

        atomicAdd(fx + ii, -FdTdRi.x);
        atomicAdd(fy + ii, -FdTdRi.y);
        atomicAdd(fz + ii, -FdTdRi.z);

        atomicAdd(fx + jj, -FdTdRj.x);
        atomicAdd(fy + jj, -FdTdRj.y);
        atomicAdd(fz + jj, -FdTdRj.z);

        atomicAdd(fx + kk, -FdTdRk.x);
        atomicAdd(fy + kk, -FdTdRk.y);
        atomicAdd(fz + kk, -FdTdRk.z);

        atomicAdd(fx + ll, -FdTdRl.x);
        atomicAdd(fy + ll, -FdTdRl.y);
        atomicAdd(fz + ll, -FdTdRl.z);
    }
    /*
    atomicAdd(sion->xx+ll, fpxx);
    atomicAdd(sion->yy+ll, fpyy);
    atomicAdd(sion->zz+ll, fpzz);
    atomicAdd(sion->xy+ll, fpxy);
    atomicAdd(sion->xz+ll, fpxz);
    atomicAdd(sion->yz+ll, fpyz);
     */
    atomicAdd(e + jj, etorstot);
}

//__global__ void torsionKernelFast(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, GPUVIRIALS *sion, int * torsionListI, int *torsionListJ, int *torsionListK, int *torsionListL, int * torsionListStarts, int * torsionListSizes, int *atmTorsionListStarts, int * atmTorsionListSizes,int *numAtomsInResiduePrefixSum, double *kchi1, double *delta1, int *n1,THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void torsionKernelFast(STATE *state, double *rx, double *ry, double *rz, gid_type *labelS, int *resIDs, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * torsionListI, int *torsionListJ, int *torsionListK, int *torsionListL, int * torsionListStarts, int * torsionListSizes, int *atmTorsionListStarts, int * atmTorsionListSizes, int *numAtomsInResiduePrefixSum, double *kchi1, double *delta1, int *n1, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int jj = r_back[pid];
    if (jj >= nLocal)
    {
        return;
    }
    long long unsigned label = labelS[pid];
    int resID = resIDs[jj];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;


    //printf("pid %i\n", pid);
    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int torsion_start_idx = atmTorsionListStarts[atm_idx];
    int torsion_end_idx = atmTorsionListSizes[atm_idx] + torsion_start_idx;

    //printf("jj %i atmidx %i  resid %i start %i end %i\n", jj,atm_idx,resID, torsion_start_idx, torsion_end_idx);
    if (torsion_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((jj >= nLocal))
    {
        return;
    }

    STATE gpustate;
    STATE * gstate = &gpustate;
    gstate->rx = rx;
    gstate->ry = ry;
    gstate->rz = rz;

    //atom J is used for determining which atom in angle owns this force term
    int indexJ = pid;

    //double fpxx=0; double fpyy=0; double fpzz=0; double fpxy=0; double fpxz=0; double fpyz=0;
    double etorstot = 0;
    //double et =0;
    for (int bi = torsion_start_idx; bi < torsion_end_idx; bi++)
    {
        //load registers and calc vectors IJ and JK 
        int indexI = pid - torsionListJ[bi] + torsionListI[bi];
        int indexK = pid - torsionListJ[bi] + torsionListK[bi];
        int indexL = pid - torsionListJ[bi] + torsionListL[bi];
        int ii = r_back[indexI];
        int kk = r_back[indexK];
        int ll = r_back[indexL];

        //Energy
        double tors;
        double sinX;

        THREE_VECTOR dTdrI;
        THREE_VECTOR dTdrJ;
        THREE_VECTOR dTdrK;
        THREE_VECTOR dTdrL;
        THREE_SMATRIX ssion;
        double kchi = kchi1[bi];
        double delta = delta1[bi];
        int n = n1[bi];

        bioDihedralFast(gstate, h1, indexI, indexJ, indexK, indexL, &tors, &sinX,
                        &dTdrI, &dTdrJ, &dTdrK, &dTdrL, &ssion);


        double etors = kchi * (1 + cos(n * tors - delta));
        etorstot = etorstot + etors;

        //printf("ii %i j %i rx1 %f rx2 %f tors %f etors %f \n", ii, r_back[indexJ],rx[indexI],rx[indexJ],  tors, etors);
        if (etors > 0.1 || etors<-0.1)
        {
            printf("I=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            //printf("Torsion =%f Torsion angle=%f kchi=%f n=%d delta=%f Torsion Energy=%f\n", 
            //    tors*RAD2DEG, tors2*RAD2DEG, torsPtr->kchi, torsPtr->n, torsPtr->delta*RAD2DEG, etors);  
        }
        //Force
        double absX = fabs(sinX);

        double ktors;

        if (absX > FLOAT_EPS)
        {
            ktors = kchi * n * sin(n * tors - delta) / sinX;
            //negative sign canceled out
        }
        else
        {
            double nX = n*tors;
            double nX2 = nX*nX;
            double nX4 = nX2*nX2;
            double nX6 = nX4*nX2;
            double nX8 = nX4*nX4;
            double nX10 = nX8*nX2;
            double X2 = tors*tors;
            double X4 = X2*X2;
            double X6 = X4*X2;
            double X8 = X4*X4;
            double X10 = X8*X2;
            double ratio = n * (1 - nX2 / 6 + nX4 / 120 - nX6 / 5040 + nX8 / 362880 - nX10 / 39916800) / (1 - X2 / 6 + X4 / 120 - X6 / 5040 + X8 / 362880 - X10 / 39916800);
            if (delta < NEAR_ZERO_ANGLE)
            { // for delta=0
                ktors = kchi * n*ratio;
            }
            else if (delta > NEAR_180_ANGLE)
            { // for delta=180
                ktors = -kchi * n*ratio;
                // sin(x-pi)=-sin(pi-x)=-sin(x)
            }
            else
            {
                //printf("Warning: torsPtr->delta doesn't equal to 0 or 180");
                ktors = kchi * n*ratio;
            }
        }


        //ktors=-ktors; // Flip the sign for bioVec
        /*
        fpxx += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.x+(FdTdRk.x-FdTdRj.x)*v_kj.x+(FdTdRk.x-FdTdRl.x)*v_kl.x);
        fpxy += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.y+(FdTdRk.x-FdTdRj.x)*v_kj.y+(FdTdRk.x-FdTdRl.x)*v_kl.y);
        fpxz += 0.5*((FdTdRi.x-FdTdRj.x)*v_ij.z+(FdTdRk.x-FdTdRj.x)*v_kj.z+(FdTdRk.x-FdTdRl.x)*v_kl.z); 
        fpyy += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.y+(FdTdRk.y-FdTdRj.y)*v_kj.y+(FdTdRk.y-FdTdRl.y)*v_kl.y);
        fpyz += 0.5*((FdTdRi.y-FdTdRj.y)*v_ij.z+(FdTdRk.y-FdTdRj.y)*v_kj.z+(FdTdRk.y-FdTdRl.y)*v_kl.z);
        fpzz += 0.5*((FdTdRi.z-FdTdRj.z)*v_ij.z+(FdTdRk.z-FdTdRj.z)*v_kj.z+(FdTdRk.z-FdTdRl.z)*v_kl.z);
         */

        THREE_VECTOR FdTdRi = bioScale(dTdrI, ktors);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, ktors);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, ktors);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, ktors);

        atomicAdd(fx + ii, -FdTdRi.x);
        atomicAdd(fy + ii, -FdTdRi.y);
        atomicAdd(fz + ii, -FdTdRi.z);

        atomicAdd(fx + jj, -FdTdRj.x);
        atomicAdd(fy + jj, -FdTdRj.y);
        atomicAdd(fz + jj, -FdTdRj.z);

        atomicAdd(fx + kk, -FdTdRk.x);
        atomicAdd(fy + kk, -FdTdRk.y);
        atomicAdd(fz + kk, -FdTdRk.z);

        atomicAdd(fx + ll, -FdTdRl.x);
        atomicAdd(fy + ll, -FdTdRl.y);
        atomicAdd(fz + ll, -FdTdRl.z);

        //atomicAdd(sion->xx+jj, -ssion.xx*ktors);
        //atomicAdd(sion->yy+jj, -ssion.yy*ktors);
        //atomicAdd(sion->zz+jj, -ssion.zz*ktors);
        //atomicAdd(sion->xy+jj, -ssion.xy*ktors);
        //atomicAdd(sion->xz+jj, -ssion.xz*ktors);
        //atomicAdd(sion->yz+jj, -ssion.yz*ktors);

        atomicAdd(sxx + jj, -ssion.xx * ktors);
        atomicAdd(syy + jj, -ssion.yy * ktors);
        atomicAdd(szz + jj, -ssion.zz * ktors);
        atomicAdd(sxy + jj, -ssion.xy * ktors);
        atomicAdd(sxz + jj, -ssion.xz * ktors);
        atomicAdd(syz + jj, -ssion.yz * ktors);
    }

    atomicAdd(e + jj, etorstot);
}





//__global__ void bpairKernel(STATE *state, double *rx, double *ry, double *rz, double *q, gid_type *labelS, int *resid, int * r_back, double *e, GPUVIRIALS *sion, int * bpairListI, int *bpairListJ, int * bpairListStarts, int * bpairListSizes, int *atmBpairListStarts, int * atmBpairListSizes,int *numAtomsInResiduePrefixSum, double *eps, double *sig, double *shifts, double rmax, double ke, double iepsilon_r, THREE_MATRIX *h1, int nLocal, int nIon) 

__global__ void bpairKernel(STATE *state, double *rx, double *ry, double *rz, double *q, gid_type *labelS, int *resid, int * r_back, double *e, double *sxx, double *syy, double *szz, double *sxy, double *sxz, double *syz, int * bpairListI, int *bpairListJ, int * bpairListStarts, int * bpairListSizes, int *atmBpairListStarts, int * atmBpairListSizes, int *numAtomsInResiduePrefixSum, double *eps, double *sig, double *shifts, double rmax, double ke, double iepsilon_r, THREE_MATRIX *h1, int nLocal, int nIon)
{
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nIon)
    {
        return;
    }
    int ii = r_back[pid];

    long long unsigned label = labelS[pid];
    int resID = resid[ii];
    int atomID = label & atmMask;
    atomID = label & atmgrpMask;

    int atm_idx = numAtomsInResiduePrefixSum[resID] + atomID;
    int bpair_start_idx = atmBpairListStarts[atm_idx];
    int bpair_end_idx = atmBpairListSizes[atm_idx] + bpair_start_idx;

    //printf("ii %i resid %i start %i\n", ii,resID, bond_start_idx);
    if (bpair_start_idx == -1)
    { //no bonds owned by atom
        return;
    }
    if ((ii >= nLocal))
    {
        return;
    }
    //printf("ii %i bond_start %i\n", ii, bond_start_idx);
    double ix = rx[pid];
    double iy = ry[pid];
    double iz = rz[pid];
    double eb = 0;

    double ffx = 0;
    double ffy = 0;
    double ffz = 0;

    for (int bi = bpair_start_idx; bi < bpair_end_idx; bi++)
    {
        double fpxx = 0;
        double fpyy = 0;
        double fpzz = 0;
        double fpxy = 0;
        double fpxz = 0;
        double fpyz = 0;
        int indexJ = pid - bpairListI[bi] + bpairListJ[bi];
        double jx = rx[indexJ];
        double jy = ry[indexJ];
        double jz = rz[indexJ];

        double bx = ix - jx;
        double by = iy - jy;
        double bz = iz - jz;

        PBC(preduceGPUO3b(h1, &bx, &by, &bz);)

            //    double kbb = kb[bi];
            double r2 = bx * bx + by * by + bz*bz;
        double rinv = rsqrt(r2);
        double r = r2*rinv;
        if (r > rmax) continue; // if the pair is larger than the cutoff, don't count it;
        // L-J term
        double sigma = sig[bi];
        double epsilon = -eps[bi];
        double shift = -shifts[bi];
        double sigma_r = sigma * rinv;
        double s2 = sigma_r*sigma_r;
        double s4 = s2*s2;
        double s6 = s4*s2;
        double s12 = s6*s6;
        double dvdr = 24.0 * epsilon * (s6 - 2.0 * s12) * rinv;
        double ebpair = 4.0 * epsilon * (s12 - s6) + shift;

        // Electrostatic term.
        double qI = q[pid];
        double qJ = q[indexJ];

        double epseduo = -ke * qI * qJ * iepsilon_r;
        double ebelec = epseduo * rinv;
        double dedr = -epseduo * rinv*rinv;

        double unit_x = bx*rinv;
        double unit_y = by*rinv;
        double unit_z = bz*rinv;

        double kforce = -(dvdr + dedr);
        double fxDelta = kforce*unit_x;
        double fyDelta = kforce*unit_y;
        double fzDelta = kforce*unit_z;

        ffx += fxDelta;
        ffy += fyDelta;
        ffz += fzDelta;
        eb += (ebelec + ebpair);

        fpxx = fxDelta*bx;
        fpxy = fxDelta*by;
        fpxz = fxDelta*bz;
        fpyy = fyDelta*by;
        fpyz = fyDelta*bz;
        fpzz = fzDelta*bz;

        int jj = r_back[indexJ];
        //atomicAdd(sion->xx+jj, -fpxx);
        //atomicAdd(sion->yy+jj, -fpyy);
        //atomicAdd(sion->zz+jj, -fpzz);
        //atomicAdd(sion->xy+jj, -fpxy);
        //atomicAdd(sion->xz+jj, -fpxz);
        //atomicAdd(sion->yz+jj, -fpyz);

        atomicAdd(sxx + jj, -fpxx);
        atomicAdd(syy + jj, -fpyy);
        atomicAdd(szz + jj, -fpzz);
        atomicAdd(sxy + jj, -fpxy);
        atomicAdd(sxz + jj, -fpxz);
        atomicAdd(syz + jj, -fpyz);

        atomicAdd(fx + jj, -fxDelta);
        atomicAdd(fy + jj, -fyDelta);
        atomicAdd(fz + jj, -fzDelta);
        //printf("ii %i jj %i r%f shift %f ep %f  elj %f eco %f fx %f fy %f fz %f\n", ii, jj,r,shift, eb,ebelec, ebpair, fxDelta, fyDelta, fzDelta);
    }
    atomicAdd(e + ii, eb);
    atomicAdd(fx + ii, ffx);
    atomicAdd(fy + ii, ffy);
    atomicAdd(fz + ii, ffz);

    return;
}

void charmmConvalentGPU(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e) 
{
    //allocate ResiCon, move to parms
    //allocResiCon(sys, parms, parms->charmmParms);	

    //grab some pointers

    GPUNLIST *gnlist=sys->collection->gnlist;
    STATE *gsh = sys->collection->gpustate_h;
    STATE *state = sys->collection->state;
    CharmmBONDGPU_PARMS * gparms = parms->gpu_bondparms_h;
    THREE_MATRIX *h1 = gnlist->hmat_g;
    CharmmANGLEGPU_PARMS * gaparms = parms->gpu_angleparms_h;
    CharmmANGLEGPU_PARMS * gcaparms = parms->gpu_cosangleparms_h;
    CharmmANGLEGPU_PARMS * grebaparms = parms->gpu_rebangleparms_h;
    CharmmTORSIONGPU_PARMS * gtparms = parms->gpu_torsionparms_h;
    CharmmIMPROPERGPU_PARMS * giparms = parms->gpu_improperparms_h;
    CharmmBPAIRGPU_PARMS * gbparms = parms->gpu_bpairparms_h;
    COLLECTION * col = sys->collection;
    if (first)
    {
        //cudaMalloc the labels
        //todo: move to cuda utils 
        gpu_allocator(gparms->labelSorted, state->nion);
        gpu_allocator(gsh->label, state->nion);

        //move labels to memory, todo we can remove the first one
        gpu_memcpy_host2device(gsh->label, state->label, state->nion);
        gpu_memcpy_host2device(gparms->labelSorted, state->label, state->nion);

        //allocate Back Pointers
        //gpu_allocator(gparms->r_backBondInp, state->nion);
        CUDA_SAFE_CALL(cudaMallocManaged((void **) &(gparms->r_backBondInp), sizeof (int)*state->nion);)
        gpu_allocator(gparms->r_backBond, state->nion);

        //fill back pointers
        for (int i = 0; i < state->nion; i++)
        {
            gparms->r_backBondInp[i] = i;
        }

        //radix sort
        radix_sort_cub(gparms->r_backBondInp, gparms->r_backBond, gsh->label, gparms->labelSorted, state->nlocal, state->nion);

        //check the results are sane

        for (int i = 0; i < state->nion; i++)
        {
            //printf("rback %i %i %lu\n", i, gparms->r_backBond[i], gparms->labelSorted[i]);
        }
        first = 0;
    }



    //printf("permute\n");
    //apply the premutation obtained from previous sort of gids
    int blockSize = SHARED_BLIST;
    int gridSize = ceil((float) state->nion / blockSize);
    permuteParticlesWithGid << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gbparms->q, gsh->label, gparms->labelSorted, gparms->r_backBond, gparms->resID, gparms->resIDs, state->nion);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    //cudaStreamSynchronize(0);

    //printf("do bond\n");

    //double *sxx = sion->xx;
    //double *syy = sion->yy;
    //double *szz = sion->zz;
    //double *sxy = sion->xy;
    //double *sxz = sion->xz;
    //double *syz = sion->yz;
    //get bond force and energies
    DEBUG_ENERGY(
                 int reduceSize = 3;
                 double* esumAll_d = NULL;
                 gpu_allocator(esumAll_d, reduceSize);
                 double *newEnergy = (double*) malloc(reduceSize * sizeof (double));
                 double oldEnergy = 0;
                 double cE = units_convert(1.0, NULL, "kJ*mol^-1");

                 size_t temp_storage_bytes = 0;
                 void *d_temp_storage = NULL;
                 size_t temp_storage_bytes2 = 0;
                 void *d_temp_storage2 = NULL;

                 int nLocal = state->nlocal;

                 cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, col->e_all, esumAll_d, nLocal);
                 cudaMalloc(&d_temp_storage, temp_storage_bytes);

                 cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, col->e_all, esumAll_d + 1, nLocal);
                 cudaMalloc(&d_temp_storage2, temp_storage_bytes2);

                 //reduce energy
                 cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, col->e_all, esumAll_d, nLocal);
                 cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, col->e_all, esumAll_d + 1, nLocal);
                 cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, col->e_all, esumAll_d + 2, nLocal);
                 gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                 printf("Before any energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                 oldEnergy = newEnergy[0];
                 )

    if (parms->bond_fcn)
    {
        //bondKernel<<<gridSize, blockSize>>>(col->gpustate, gparms->rxs,gparms->rys,gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all,gnlist->gpu_sion,  gparms->bondListI,  gparms->bondListJ, gparms->resiBondListStarts, gparms->resiBondListSizes, gparms->atmBondListStarts, gparms->atmBondListSizes, gparms->numAtomsInResiduePrefixSum, gparms->b0, gparms->kb,h1, state->nlocal, state->nion);
        bondKernel << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gparms->bondListI, gparms->bondListJ, gparms->resiBondListStarts, gparms->resiBondListSizes, gparms->atmBondListStarts, gparms->atmBondListSizes, gparms->numAtomsInResiduePrefixSum, gparms->b0, gparms->kb, h1, state->nlocal, state->nion);

        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Bond energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )

    }
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
        //cudaStreamSynchronize(0);


    if (parms->angle_fcn)
    {
        //angleKernel<<<gridSize, blockSize>>>(col->gpustate,gparms->rxs,gparms->rys,gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all,gnlist->gpu_sion, gaparms->angleListI, gaparms->angleListJ, gaparms->angleListK, gaparms->resiAngleListStarts, gaparms->resiAngleListSizes, gaparms->atmAngleListStarts, gaparms->atmAngleListSizes,gparms->numAtomsInResiduePrefixSum, gaparms->ktheta, gaparms->theta0, gaparms->kub, gaparms->s0,h1, state->nlocal, state->nion); 
        angleKernel << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gaparms->angleListI, gaparms->angleListJ, gaparms->angleListK, gaparms->resiAngleListStarts, gaparms->resiAngleListSizes, gaparms->atmAngleListStarts, gaparms->atmAngleListSizes, gparms->numAtomsInResiduePrefixSum, gaparms->ktheta, gaparms->theta0, gaparms->kub, gaparms->s0, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Angle energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )
    }


    if (parms->cosangle_fcn)
    {
        //angleCosineKernel<<<gridSize, blockSize>>>(col->gpustate,  gparms->rxs,gparms->rys,gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, col->e_all, col->gpu_sion,gcaparms->angleListI, gcaparms->angleListJ, gcaparms->angleListK, gcaparms->resiAngleListStarts, gcaparms->resiAngleListSizes, gcaparms->atmAngleListStarts, gcaparms->atmAngleListSizes,gparms->numAtomsInResiduePrefixSum, gcaparms->ktheta, gcaparms->theta0, gcaparms->kub, gcaparms->s0, h1, state->nlocal, state->nion); 
        angleCosineKernel << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gcaparms->angleListI, gcaparms->angleListJ, gcaparms->angleListK, gcaparms->resiAngleListStarts, gcaparms->resiAngleListSizes, gcaparms->atmAngleListStarts, gcaparms->atmAngleListSizes, gparms->numAtomsInResiduePrefixSum, gcaparms->ktheta, gcaparms->theta0, gcaparms->kub, gcaparms->s0, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Cosine Angle energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )
            CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }

    if (parms->rebangle_fcn)
    {
        angleRebKernel << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, grebaparms->angleListI, grebaparms->angleListJ, grebaparms->angleListK, grebaparms->resiAngleListStarts, grebaparms->resiAngleListSizes, grebaparms->atmAngleListStarts, grebaparms->atmAngleListSizes, gparms->numAtomsInResiduePrefixSum, grebaparms->ktheta, grebaparms->theta0, grebaparms->kub, grebaparms->s0, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Restricted Bending Angle energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )
            CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }


    if (parms->torsion_fcn)
    {
        //torsionKernelFast<<<gridSize, blockSize>>>(col->gpustate,  gparms->rxs,gparms->rys,gparms->rzs, gparms->labelSorted,  gparms->resID, gparms->r_backBond, col->e_all,col->gpu_sion, gtparms->torsionListI, gtparms->torsionListJ, gtparms->torsionListK, gtparms->torsionListL, gtparms->resiTorsionListStarts, gtparms->resiTorsionListSizes, gtparms->atmTorsionListStarts, gtparms->atmTorsionListSizes,gparms->numAtomsInResiduePrefixSum, gtparms->kchi, gtparms->delta, gtparms->n, h1, state->nlocal, state->nion);
        torsionKernelFast << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gtparms->torsionListI, gtparms->torsionListJ, gtparms->torsionListK, gtparms->torsionListL, gtparms->resiTorsionListStarts, gtparms->resiTorsionListSizes, gtparms->atmTorsionListStarts, gtparms->atmTorsionListSizes, gparms->numAtomsInResiduePrefixSum, gtparms->kchi, gtparms->delta, gtparms->n, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Torsion energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )
            CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }



    if (parms->improper_fcn)
    {
        //improperKernelFast<<<gridSize, blockSize>>>(col->gpustate,  gparms->rxs,gparms->rys,gparms->rzs, gparms->labelSorted,  gparms->resID, gparms->r_backBond, col->e_all,col->gpu_sion, giparms->improperListI, giparms->improperListJ, giparms->improperListK, giparms->improperListL, giparms->resiImproperListStarts, giparms->resiImproperListSizes, giparms->atmImproperListStarts, giparms->atmImproperListSizes,gparms->numAtomsInResiduePrefixSum, giparms->kpsi, giparms->psi0, h1, state->nlocal, state->nion);
        improperKernelFast << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, giparms->improperListI, giparms->improperListJ, giparms->improperListK, giparms->improperListL, giparms->resiImproperListStarts, giparms->resiImproperListSizes, giparms->atmImproperListStarts, giparms->atmImproperListSizes, gparms->numAtomsInResiduePrefixSum, giparms->kpsi, giparms->psi0, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Improper energy: %f     min %f   max %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2]);
                     oldEnergy = newEnergy[0];
                     )
            CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }

    if (parms->bpair_fcn)
    {

        double ke = parms->gpu_ljparms_h->ke;
        double iepsilon_r = parms->gpu_ljparms_h->iepsilon_r;
        //bpairKernel<<<gridSize, blockSize>>>(col->gpustate,  gparms->rxs,gparms->rys,gparms->rzs,gbparms->q, gparms->labelSorted, gparms->resID, gparms->r_backBond, col->e_all, col->gpu_sion, gbparms->bpairListI, gbparms->bpairListJ, gbparms->resiBpairListStarts, gbparms->resiBpairListSizes, gbparms->atmBpairListStarts, gbparms->atmBpairListSizes,gparms->numAtomsInResiduePrefixSum, gbparms->eps, gbparms->sigma, gbparms->shift, parms->rmax, ke, iepsilon_r, h1, state->nlocal, state->nion); 
        bpairKernel << <gridSize, blockSize>>>(col->gpustate, gparms->rxs, gparms->rys, gparms->rzs, gbparms->q, gparms->labelSorted, gparms->resID, gparms->r_backBond, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gbparms->bpairListI, gbparms->bpairListJ, gbparms->resiBpairListStarts, gbparms->resiBpairListSizes, gbparms->atmBpairListStarts, gbparms->atmBpairListSizes, gparms->numAtomsInResiduePrefixSum, gbparms->eps, gbparms->sigma, gbparms->shift, parms->rmax, ke, iepsilon_r, h1, state->nlocal, state->nion);
        DEBUG_ENERGY(
                     //reduce energy
                     cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
                     cub::DeviceReduce::Min(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 1, nLocal);
                     cub::DeviceReduce::Max(d_temp_storage2, temp_storage_bytes2, gnlist->e_all, esumAll_d + 2, nLocal);
                     gpu_memcpy_device2host(newEnergy, esumAll_d, 3);
                     printf("Pair energy: %f     min %f   max %f   cE %f\n", cE * (newEnergy[0] - oldEnergy), newEnergy[1], newEnergy[2], cE);
                     )
            CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }

    DEBUG_ENERGY(
                 if (d_temp_storage) cudaFree(d_temp_storage);
                 if (d_temp_storage2) cudaFree(d_temp_storage2);
                 if (esumAll_d) cudaFree(esumAll_d);
                 if (newEnergy) free(newEnergy);
                 )
                    //temporary: send energies to cpu and sum
                    /*
                    double * e_cpu = (double *)malloc(sizeof(double)*state->nlocal);
                    gpu_memcpy_device2host(e_cpu,  col->e_all, state->nlocal);
                    double e_sum =0;

                    for (int i=0 ; i<state->nlocal;i++){
                       //printf("i %i e_cpu %f \n", i, e_cpu[i]);
                       e_sum+=e_cpu[i];
                         }
                    //e->eion+=e_sum;
                    //printf("bond esum %f \n", e_sum);
                     */
                }



#define radix_groupSize 16 //threads in a radix group
#define radix_bitL 8 //bits used in radi sort
#define radix_blocks 16 //usually equal to number of MPs
#define radix_num 256 //2^radix_bitL 
#define radix_counterNum //num counters = radix_groups*blockSize*numradices
#define radix_groups 12 //determined empirically
#define radix_blockSize 192 //radix block size

__global__ void countRadices(gid_type* gids, int *ids, int *counts, int pass, int nLocal, int nIon)
{
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int groupId = tid / radix_groupSize;
    if (pid >= nIon) return; //guard statement to avoid segfaults

    //we'd ideally like a  set of counters for each thread, but
    //that would require radix_num*blockSize*numblocks of shared memory - too much
    //so instead, we divide each blocks into 'groups', such that each group
    //now has a set of counters
    __shared__ int radix_counters[radix_num * radix_groups];
    for (int i = pid; i < radix_num * radix_groups; i += blockDim.x)
    {
        radix_counters[i] = 0;
    }
    __syncthreads();

    //each block gets assigned a portion of the gid array
    int portion_size = max(nIon / radix_blocks, nIon);
    int start = portion_size * blockIdx.x;
    int end = min(portion_size * (blockIdx.x + 1), nIon);
    int mask = 1; //(1 << radix_bitL)-1;
    mask = mask << radix_bitL*pass; //shift the mask to analyze bits for current pass
    printf("pid %i gid %i tid %i s/e %i/%i\n", pid, groupId, tid, start, end);
    for (int i = start; i < end; i += radix_groupSize)
    {
        //combine appropriate mask to check if radix exists, atomic add to count 
        //occurences
        int result = gids[i] & mask;
        printf("pid %i mask %i gid % result %i\n", pid, mask, gids[i], result);
        atomicAdd(radix_counters + groupId * radix_groupSize + result, 1);
    }
    __syncthreads();

    //commit to global memory
    if (pid < radix_num)
    {
        for (int i = pid; i < radix_num; i += nIon)
        {
            for (int j = 0; j < radix_groups; j++)
            {
                atomicAdd(counts + radix_groupSize + j, radix_counters[j]);
            }
        }
    }
}
#if 0
int *d_key_buf; // e.g., [8, 6, 7, 5, 3, 0, 9]
int *d_key_alt_buf; // e.g., [        ...        ]
int *d_value_buf; // e.g., [0, 1, 2, 3, 4, 5, 6]
int *d_value_alt_buf; // e.g., [        ...        ]


int n = 4096;
int num_items = n;
int *h_key_buf = (int*) malloc(sizeof (int)*n);
int *h_value_buf = (int*) malloc(sizeof (int)*n);
int *h_key_buf2 = (int*) malloc(sizeof (int)*n);
int *h_value_buf2 = (int*) malloc(sizeof (int)*n);


for (int i = 0; i < n; i++)
{
    h_key_buf[i] = n - i - 1;
    h_value_buf[i] = i;
}

gpu_allocator(d_key_buf, n);
gpu_allocator(d_value_buf, n);
gpu_allocator(d_key_alt_buf, n);
gpu_allocator(d_value_alt_buf, n);

printf("done allocating\n");
gpu_memcpy_host2device(d_key_buf, h_key_buf, n);
gpu_memcpy_host2device(d_value_buf, h_value_buf, n);

printf("done memcpying\n");
// Create a set of DoubleBuffers to wrap pairs of device pointers
cub::DoubleBuffer<int> d_keys(d_key_buf, d_key_alt_buf);
cub::DoubleBuffer<int> d_values(d_value_buf, d_value_alt_buf);
printf("pointers wrapped\n");
// Determine temporary device storage requirements
void *d_temp_storage = NULL;
size_t temp_storage_bytes = 0;
cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, num_items);
// Allocate temporary storage
printf("temp storage\n");
cudaMalloc(&d_temp_storage, temp_storage_bytes);
// Run sorting operation
cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, num_items);
// d_keys.Current()      <-- [9, 8, 7, 6, 5, 3, 0]
gpu_memcpy_device2host(h_key_buf2, d_key_buf, n);
gpu_memcpy_device2host(h_value_buf2, d_value_buf, n);

for (int i = 0; i < n; i++)
{
    //printf("i %i key %i value %i\n",i, h_key_buf[i], h_value_buf2[i]);

#endif

    //int * bondListIg,* bondListJg;
    //cudaMalloc((void **) &bondListIg, sizeof(int)*bondListSize); 
    //cudaMalloc((void **) &bondListJg, sizeof(int)*bondListSize); 
    //cudaMemPrefetchAsync(bondListIg, bondListI, sizeof(int)*bondListSize, 
    //cudaMemcpy(bondListIg, bondListI, sizeof(int)*bondListSize, cudaMemcpyHostToDevice);
    //cudaMemcpy(bondListJg, bondListJ, sizeof(int)*bondListSize, cudaMemcpyHostToDevice);



    /*
       gpu_allocator(gpu_parms->angleListI, bondListSize);
       gpu_allocator(gpu_parms->angleListJ, bondListSize);
       gpu_allocator(gpu_parms->angleListK, bondListSize);

       gpu_allocator(gpu_parms->torsListI, bondListSize);
       gpu_allocator(gpu_parms->torsListJ, bondListSize);
       gpu_allocator(gpu_parms->torsListK, bondListSize);
       gpu_allocator(gpu_parms->torsListK, bondListSize);
     */
    /*
    int *bondListI  = (int *)malloc(sizeof(int)*bondListSize);
         int *bondListJ  = (int *)malloc(sizeof(int)*bondListSize);
    int *angleListI = (int *)malloc(sizeof(int)*angleListSize);
         int *angleListJ = (int *)malloc(sizeof(int)*angleListSize);
         int *angleListK = (int *)malloc(sizeof(int)*angleListSize);
         int *torsListI  = (int *)malloc(sizeof(int)*torsListSize);
         int *torsListJ  = (int *)malloc(sizeof(int)*torsListSize);
         int *torsListK  = (int *)malloc(sizeof(int)*torsListSize);
         int *torsListL  = (int *)malloc(sizeof(int)*torsListSize);
     */

    /*
    int * bondListIg, bondListJg;
    gpu_allocator(bondListIg, bondListSize);
    gpu_allocator(bondListJg, bondListSize);
    gpu_memcpy_host2device(bondListIg,  bondListI, bondListSize);
     */

    //allocate memory that will give us offsets into a bond/angle/torsion list
    //for a particular residue (within bond/angle/torsion list buffers)
    /*
         int * bondListStarts  = (int *) malloc(sizeof(charmmParms->resiConnSize));
         int * bondListSizes    = (int *) malloc(sizeof(charmmParms->resiConnSize));
         int * angleListStarts = (int *) malloc(sizeof(charmmParms->resiConnSize));
         int * angleListSizes   = (int *) malloc(sizeof(charmmParms->resiConnSize));
         int * torsListStarts  = (int *) malloc(sizeof(charmmParms->resiConnSize));
         int * torsListSizes    = (int *) malloc(sizeof(charmmParms->resiConnSize));
     */


    //printf("ii %i jj %i ee %f fx %f fy %f fz %f rij %f costheta %f cos %f\n", ii, jj, ee, fxDeltaK, fyDeltaK, fzDeltaK,ijr, kx, cosTheta);
    /*     
              atomicAdd(fx+ii, fxDeltaI);
              atomicAdd(fy+ii, fyDeltaI);
              atomicAdd(fz+ii, fzDeltaI);

              atomicAdd(fx+kk, fxDeltaK);
              atomicAdd(fy+kk, fyDeltaK);
              atomicAdd(fz+kk, fzDeltaK);

              fxj-= (fxDeltaI + fxDeltaK);
              fyj-= (fyDeltaI + fyDeltaK);
              fzj-= (fzDeltaI + fzDeltaK);
     */

    // Pressure
    /*
              double  fpxx = (fxDeltaI*vec_ij.x+fxDeltaK*vec_kj.x);
              double  fpxy = (fxDeltaI*vec_ij.y+fxDeltaK*vec_kj.y);
              double  fpxz = (fxDeltaI*vec_ij.z+fxDeltaK*vec_kj.z);
              double  fpyy = (fyDeltaI*vec_ij.y+fyDeltaK*vec_kj.y);
              double  fpyz = (fyDeltaI*vec_ij.z+fyDeltaK*vec_kj.z);
              double  fpzz = (fzDeltaI*vec_ij.z+fzDeltaK*vec_kj.z);   
     */


    /*
         THREE_VECTOR *v_ij = &v_ij1;
         THREE_VECTOR *v_kj = &v_kj1;
         THREE_VECTOR *v_kl = &v_kl1;
         bioVec(gstate, indexI, indexJ, v_ij);
         bioVec(state, indexK, indexJ, v_kj);
         bioVec(state, indexK, indexL, v_kl);

         THREE_VECTOR v_A=cross_g(v_ij, v_kj);
         THREE_VECTOR v_B=cross_g(v_kj, v_kl);       
         THREE_VECTOR v_C=cross_g(&v_A, v_kj);    
        
         double v1_A=bioNorm(v_A);
         double v1_B=bioNorm(v_B);
         double v1_C=bioNorm(v_C);
    
         if(v1_A<1.0e-7){
            printf("Warning: A vector rij X rkj is too small\n");
         }
         if(v1_B<1.0e-7){
            printf("Warning: B vector rkj X rkl is too small\n");
         }    
         if(v1_C<1.0e-7){
            printf("Warning: C vector rkj X A is too small\n");
         }    
    
    
         double dotCB=dot1_g(v_C, v_B);
         double sinX=dotCB/(v1_C*v1_B); 
     
           
         double dotAB=dot1_g(v_A, v_B);
         double cosX=dotAB/(v1_A*v1_B);
     
         double angleX;     

       if(cosX>0){
            angleX=atan(sinX/cosX);
        }else if(sinX>=0 && cosX<0){
            angleX=atan(sinX/cosX)+M_PI;
        }else if(sinX<0  && cosX<0){
            angleX=atan(sinX/cosX)-M_PI;
        }else if(sinX>0  && cosX==0){
            angleX=M_PI/2;
        }else if(sinX<0  && cosX==0){
            angleX=-M_PI/2;
        }else if(sinX==0  && cosX==0){
            //normally won't reach here.
            printf("ERROR: undefined angle\n");
            angleX=0;
         }

        double diffsinX=fabs(sinX-sin(angleX));
        if(diffsinX>1e-4){
            printf("ERROR: X angle %f is not consist with sinX, %f \n", angleX, sinX);
        }   

        double coef=1/(v1_A*v1_B);
        double coef_a=-cosX/(v1_A*v1_A);
        double coef_b=-cosX/(v1_B*v1_B);    
        //      1   B     A               1  
        // a = ---(--- - --- cosX) * (- ----)
        //     |A| |B|   |A|            sinX
        // b are similar. (-1/sinX is from dX/dcosX)
        THREE_VECTOR v_a=bioVecAdd(coef, v_B, coef_a, v_A);
        THREE_VECTOR v_b=bioVecAdd(coef, v_A, coef_b, v_B);    
    
        THREE_VECTOR dcosXdrI;
        THREE_VECTOR dcosXdrJ;
        THREE_VECTOR dcosXdrK;
        THREE_VECTOR dcosXdrL;
        // dcosX/drI = - r_kj X a
        dcosXdrI=cross_g(&v_a, v_kj);
    
        // dcosX/drL = - r_kj X b   
        dcosXdrL=cross_g(&v_b, v_kj);
    
        // dcosX/drJ = - r_ik X a + r_kl X b
        THREE_VECTOR v_ik;
        bioVec(state, indexI, indexK, &v_ik);    
        THREE_VECTOR va_ik=cross_g(&v_ik, &v_a);
        THREE_VECTOR vb_kl=cross_g(v_kl, &v_b);
        (dcosXdrJ)=bioVecAdd(-1, va_ik, 1, vb_kl);
    
        // dcosX/drK = - r_jl X b + r_ij X a
        THREE_VECTOR v_jl;
        bioVec(state, indexJ, indexL, &v_jl);    
        THREE_VECTOR vb_jl=cross_g(&v_jl, &v_b);
        THREE_VECTOR va_ij=cross_g(v_ij, &v_a);
        dcosXdrK=bioVecAdd(-1, vb_jl, 1, va_ij); 
     */




    /*
         double ix = rx[indexI];
         double iy = ry[indexI];
         double iz = rz[indexI];
     
         double ijx = ix-jx;
         double ijy = iy-jy;
         double ijz = iz-jz;

         double kx = rx[indexK];
         double ky = ry[indexK];
         double kz = rz[indexK];
 
         double lx = rx[indexL];
         double ly = ry[indexL];
         double lz = rz[indexL];    

         double kjx = kx-jx;
         double kjy = ky-jy;
         double kjz = kz-jz;

              //calculate distance between particles I,J and J,K
         double ijr2 = ijx*ijx+ijy*ijy+ijz*ijz;
         double rinv_ij = 1/sqrt(ijr2);// rsqrt(ijr2);  
              double ijr = sqrt(ijr2);//ijr2*rinv_ij;

         double kjr2 = kjx*kjx+kjy*kjy+kjz*kjz;
              double rinv_kj = 1/sqrt(kjr2);//rsqrt(kjr2); 
         double kjr =sqrt(kjr2);// kjr2*rinv_kj;

         //normalize IJ and JK vectors to get unit vectors
              ijx = ijx*rinv_ij;
              ijy = ijy*rinv_ij;
              ijz = ijz*rinv_ij;

              kjx = kjx*rinv_kj;
              kjy = kjy*rinv_kj;
              kjz = kjz*rinv_kj;

     */



