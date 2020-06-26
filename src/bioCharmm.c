#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <limits.h>
#include <mpi.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "codata.h"
#include "units.h"
#include "mpiUtils.h"
#include "ddcMalloc.h"
#include "bioCharmm.h"
#include "bioCharmmTop.h"
#include "bioCharmmPar.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "bioGid.h"
#include "ewald.h"
#include "utilities.h"
#include "system.h"
#include "ddc.h"
#include "three_algebra.h"
#include "preduce.h"
#include "accelerator.h"
#include "simulate.h"
#include "HAVEGPU.h"

/* Single precision accuracy */
#define FLOAT_EPS    1e-08

#define NEAR_ZERO_ANGLE 0.017453292519943295 // 1 degree
#define NEAR_180_ANGLE  3.12413936106985     // 179 degree


// Unit conversion
// 1.0 ry (Rydberg) = 
//    = 0.5 hartree
//    = 1.0 ry
//    = 13.605691930242388 eV
//    = 1312.7496997450642 kJ/mol
//    = 313.75470835207074 kcal/mol
//    = 2.17987197e-21 kJ
//    = 5.21001904875717e-22 kcal
//    = 2.1798719700000002e-11 erg

void printBioEnergies(BIOENERGIES en)
{
    double cE = units_convert(1.0, NULL, "kJ*mol^-1");
    double covalentE = en.bond + en.angle + en.ub + en.torsion + en.impr + en.cmap;
    double nonBondE = en.nonBond + (en.blj + en.bele);
    double totE = covalentE + nonBondE;
    printf("Covalent nonBond total (kJ/mole)\n");
    printf("%12.6f %12.6f %15.7f\n", covalentE*cE, nonBondE*cE, totE * cE);
    /*   

       //check if we want to grab gpu parms
       int use_gpu;
       object_get((OBJECT *) potential, "use_gpu", &use_gpu, INT, 1, "0");    
       if (use_gpu)
       {
          if (getRank(0) == 0) printf("using gpu martini parms\n");
          martiniNonBondGPUParms(parms);
          potential->use_gpu_list = 1;
          potential->eval_potential = martiniGPU1;
       }
       else {
          if (getRank(0) == 0) printf("using cpu martini parms\n");
       }

       timestamp("END   Bio Martini Setup");
       return parms;

       return parms;
     */
    printf("    %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n", "BOND", "ANGLE", "UrayBradly", "TORSION", "IMPROPER", "CMAP", "LJ", "ELE");
    printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f", en.bond*cE, en.angle*cE, en.ub*cE, en.torsion*cE, en.impr*cE, en.cmap * cE);
    printf(" %12.6f %12.6f\n", (en.lj + en.blj) * cE, (en.ele + en.bele) * cE);
    //printf(" %12.6f %12.6f\n", (en.lj)*cE, (en.ele)*cE);     
    //printf(" %12.6f %12.6f\n", (-en.blj)*cE, (-en.bele)*cE);     
}

typedef struct listab_struct
{
    int a, b;
} LISTAB;

typedef struct indexLabel_st
{
    int index;
    gid_type label;
} INDEXLABEL;

int compareIndexLabel(const void*pA, const void*pB)
{
    INDEXLABEL *p1 = (INDEXLABEL *) pA;
    INDEXLABEL *p2 = (INDEXLABEL *) pB;
    if (p1->label > p2->label) return 1;
    if (p1->label < p2->label) return -1;
    return 0;
}

double CHLennardJones(CharmmLJ_PARMS *ljparms, double r, double *dvdr)
{
    // CHARMM LJ: V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
    double sigma = ljparms->sigma;
    double eps = ljparms->eps;
    double ir = 1 / r;
    double sigma_r = sigma * ir;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;
    *dvdr = 12.0 * eps * (s6 - s12) * ir;
    double e = eps * (s12 - 2.0 * s6) + ljparms->shift;
    return e;
}

double CHLennardJones_setShift(double sigma, double eps, double rcut)
{
    double sigma_r = sigma / rcut;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;
    double shift = -1.0 * eps * (s12 - 2.0 * s6);
    return shift;
}
/*
CharmmLJGPU_PARMS *CharmmLJGPU_parms(CHARMMPOT_PARMS *parms)
{
    CharmmLJ_PARMS **ljParms = (CharmmLJ_PARMS **) ddcMalloc(sizeof (CharmmLJ_PARMS*) * nspecies * nspecies);
    object_get((OBJECT*) potential, "cutoff", &cutoff, WITH_UNITS, 1, "12.0", "Angstrom",NULL);
    parms->fcn = (double(*)(void *, double, double *))CHLennardJones;
    parms->rmax = cutoff;


}
 */

/*
CharmmLJ_PARMS **CharmmLJ_parms(CHARMMPOT_PARMS *parms) {
    POTENTIAL *potential = parms->potential;
    int nspecies = parms->nspecies; // this is actually nspecies=charmmParms->charmmPar->ljParmSize;
    double cutoff;
    CharmmLJ_PARMS **ljParms = (CharmmLJ_PARMS **) ddcMalloc(sizeof (CharmmLJ_PARMS*) * nspecies * nspecies);
    object_get((OBJECT*) potential, "cutoff", &cutoff, WITH_UNITS, 1, "12.0", "Angstrom",NULL);
    parms->fcn = (double(*)(void *, double, double *))CHLennardJones;
    parms->rmax = cutoff;


}
 */
CharmmLJ_PARMS CharmmLJ_parms(POTENTIAL *potential, LJCH_PARMS *speciesA, LJCH_PARMS *speciesB, double cutoff)
{
    int lkey = strlen(speciesA->atmTypeI) + strlen(speciesB->atmTypeI);

    char keyword[lkey + 16];
    double eps, sigma;
    sprintf(keyword, "eps_%s-%s", speciesA->atmTypeI, speciesB->atmTypeI);
    int flag = object_testforkeyword((OBJECT*) potential, keyword);
    if (flag == 0) sprintf(keyword, "eps_%s-%s", speciesB->atmTypeI, speciesA->atmTypeI);

    object_get_with_units((OBJECT*) potential, keyword, &eps, WITH_UNITS, 1, "-1.0", "e", NULL);
    if (eps < 0.0) eps = sqrt(speciesA->epsilon * speciesB->epsilon);

    sprintf(keyword, "sigma_%s-%s", speciesA->atmTypeI, speciesB->atmTypeI);
    flag = object_testforkeyword((OBJECT*) potential, keyword);
    if (flag == 0) sprintf(keyword, "sigma_%s-%s", speciesB->atmTypeI, speciesA->atmTypeI);

    object_get_with_units((OBJECT*) potential, keyword, &sigma, WITH_UNITS, 1, "-1.0", "l", NULL);
    if (sigma < 0.0) sigma = speciesA->rmin + speciesB->rmin;

    CharmmLJ_PARMS ljParms_ab;
    ljParms_ab.eps = eps;
    ljParms_ab.sigma = sigma;
    ljParms_ab.rcut = cutoff;
    ljParms_ab.shift = CHLennardJones_setShift(sigma, eps, cutoff);
    //ljParms_ab.rcut = sigma*cutoff;
    return ljParms_ab;
}

/*
 * Converts array of structs of LJ parms to a struct of arrays to send to gpu
 * Then sends this info to the gpu
 */
void charmmLJGPUParms(CHARMMPOT_PARMS *parms)
{
    printf("gpu parms\n\n");
    CharmmLJGPU_PARMS* ljgpu_parms = ddcMalloc(sizeof (CharmmLJGPU_PARMS));
    parms->gpu_ljparms_h = ljgpu_parms;

    SYSTEM *sys = NULL;
    sys = system_getSystem(sys);
    int nspecies = sys->nspecies;
    int nspecies2 = nspecies*nspecies;
    double *rcut = ddcMalloc(sizeof (double)*nspecies2);
    double *eps = ddcMalloc(sizeof (double)*nspecies2);
    double *sigma = ddcMalloc(sizeof (double)*nspecies2);

    //LJCH_PARMS** ljchPtr=parms->charmmParms->charmmPar->ljParms;
    //int nlist = (nspecies * nspecies + nspecies) / 2;
    CharmmLJ_PARMS ** ljparms = (CharmmLJ_PARMS **) parms->parms;
    for (int i = 0; i < sys->nspecies; i++)
    {
        eps[i] = ljparms[i]->eps;
        sigma[i] = ljparms[i]->sigma;
        rcut[i] = ljparms[i]->rcut;
    }

//allocate gpu parms struct on gpu
GPUCODE(
    cudaMalloc((void **) &(parms->gpu_ljparms), sizeof (CharmmLJGPU_PARMS));

    //allocate struct member arrays on gpu
    cudaMalloc((void **) &(ljgpu_parms->eps), sizeof (double)*nspecies);
    cudaMalloc((void **) &(ljgpu_parms->sigma), sizeof (double)*nspecies);
    cudaMalloc((void **) &(ljgpu_parms->rcut), sizeof (double)*nspecies);
    printf("\n\nmallocing new array \n");
    CUDA_SAFE_CALL(cudaMalloc((void **) &(ljgpu_parms->cg_species_index), sizeof (int)*sys->nion);)
    CUDA_SAFE_CALL(cudaMalloc((void **) &(ljgpu_parms->cg_species_index_b), sizeof (int)*sys->nion);)

    //memcpy struct members to gpu
    cudaMemcpy(ljgpu_parms->eps, eps, sizeof (double)*nspecies, cudaMemcpyHostToDevice);
    cudaMemcpy(ljgpu_parms->sigma, sigma, sizeof (double)*nspecies, cudaMemcpyHostToDevice);
    cudaMemcpy(ljgpu_parms->rcut, rcut, sizeof (double)*nspecies, cudaMemcpyHostToDevice);

    //memcpy struct to gpu
    cudaMemcpy(parms->gpu_ljparms, ljgpu_parms, sizeof (CharmmLJGPU_PARMS), cudaMemcpyHostToDevice);
)

    return;
}

CHARMMPOT_PARMS *charmm_parms(POTENTIAL *potential)
{

    timestamp("Start Bio CHARMM Setup");
    CHARMMPOT_PARMS* parms = ddcMalloc(sizeof (CHARMMPOT_PARMS));
    parms->charmmParms = ddcMalloc(sizeof (CHARMM_PARMS));
    char *topfile;
    char *parfile;
    object_get((OBJECT *) potential, "topfile", &topfile, STRING, 1, "top_all22_prot.inp");
    object_get((OBJECT *) potential, "parfile", &parfile, STRING, 1, "par_all22_prot.inp");

    if (parseCharmmParms(topfile, parfile, parms->charmmParms) < 0)
    {
        printf("Cannot load CHARMM topology or parameter files on task %d in bioCharmm\n", getRank(0));
        exit(3);
    }
    int nspecies = parms->nspecies = parms->charmmParms->charmmPar->ljParmSize;
    parms->potential = potential;
    double cutoff;

    object_get((OBJECT*) potential, "cutoff", &cutoff, WITH_UNITS, 1, "12.0", "Angstrom", NULL);
    parms->nonBond_fcn = (double(*)(void *, double, double *))CHLennardJones;
    parms->rmax = cutoff;

    LISTAB list[(nspecies * nspecies + nspecies) / 2];
    int nlist = 0;
    for (int a = 0; a < nspecies; a++)
    {
        list[nlist].a = a;
        list[nlist].b = a;
        nlist++;
    }
    for (int a = 0; a < nspecies; a++)
    {
        for (int b = 0; b < a; b++)
        {

            list[nlist].a = a;
            list[nlist].b = b;
            nlist++;
        }
    }
    CharmmLJ_PARMS *buffer = ddcMalloc(nlist * sizeof (CharmmLJ_PARMS));
    CharmmLJ_PARMS **ljParms = (CharmmLJ_PARMS **) ddcMalloc(sizeof (CharmmLJ_PARMS*) * nspecies * nspecies);
    LJCH_PARMS** ljchPtr = parms->charmmParms->charmmPar->ljParms;
    for (int k = 0; k < nlist; k++)
    {
        int a = list[k].a;
        int b = list[k].b;
        int ab = a + b*nspecies;
        int ba = b + a*nspecies;
        LJCH_PARMS *speciesA = ljchPtr[a];
        LJCH_PARMS *speciesB = ljchPtr[b];
        buffer[k] = CharmmLJ_parms(potential, speciesA, speciesB, cutoff);
        ljParms[ab] = ljParms[ba] = buffer + k;
        //if (parms->rmax < ljParms[ab]->rcut) parms->rmax = ljParms[ab]->rcut;
    }
    parms->parms = (void*) ljParms;


    object_get((OBJECT*) potential, "rmax4all", &(parms->rmax4all), DOUBLE, 1, "6.0");
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

    SYSTEM *sys = NULL;
    sys = system_getSystem(sys);
    //look for species that determine task->Residue ownership
    int isResidueOwner = 0;
    int *isResidueOwnerArr = (int*) ddcMalloc(sizeof (int)*sys->nspecies);
    object_get((OBJECT *) potential, "useAutoResidueOwnership", &(parms->useAutoResidueOwnership), INT, 1, "1");
    for (int i = 0; i < sys->nspecies; i++)
    {
        char * spname = sys->species[i]->name;
        OBJECT * spTemp = object_find(spname, "SPECIES");
        object_get((OBJECT *) spTemp, "isResidueOwner", &isResidueOwner, INT, 1, "0");
        isResidueOwnerArr[i] = isResidueOwner;
    }
    parms->isResidueOwnerSpecies = isResidueOwnerArr;

    //Allow molecules to transform?

    int transform = 0;
    object_get((OBJECT *) potential, "use_transform", &transform, INT, 1, "0");
    parms->weightedCharmm = transform;
    parms->charmmParms->charmmWeights = NULL;
    if (transform)
    {
        parms->charmmParms->charmmWeights = (BIOWEIGHTS*) ddcMalloc(4 * sizeof (BIOWEIGHTS));
        charmmTransformInit(parms->charmmParms->charmmWeights);
    }

    // Set up function pointers for connectivity energy terms
    if (transform)
    {

        parms->bond_fcn = (double (*)(void *, void *, int, int, void *, void *, void *, void *))resBondSortedWeighted;
        parms->angle_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resAngleSortedWeighted;
        parms->cosangle_fcn = NULL;
        parms->rebangle_fcn = NULL;
        parms->ureybradley_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resUreyBradleySortedWeighted;
        parms->torsion_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resTorsionSortedWeighted;
        parms->improper_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resImproperSortedWeighted;
        //      parms->bpair_fcn=resBpairSortedWeighted;
        parms->cmap_fcn = (double (*)(void *, void *, void *, int, int, void *, void *, void *))resCmap;

    }
    else
    {
        parms->bond_fcn = (double (*)(void *, void *, int, int, void *, void *, void *, void *))resBondSorted;
        parms->angle_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resAngleSorted;
        parms->cosangle_fcn = NULL;
        parms->rebangle_fcn = NULL;
        parms->ureybradley_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resUreyBradleySorted;
        parms->torsion_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resTorsionSorted;
        parms->improper_fcn = (double (*)(void *, void *, int, void *, void *, void *, void *))resImproperSorted;
        parms->cmap_fcn = (double (*)(void *, void *, void *, int, int, void *, void *, void *))resCmap;

    }


    //check if we want to grab gpu parms
    ACCELERATOR *accelerator = NULL;
    accelerator = accelerator_getAccelerator(accelerator);
    if (accelerator->itype == GPU_CUDA)
    {
        printf("use gpu\n");
        charmmLJGPUParms(parms);
    }

    timestamp("END   Bio CHARMM Setup");

    return parms;
}

RCUT_TYPE *charmmCutoff(SYSTEM *sys, CHARMMPOT_PARMS *parms, int *n)
{
    // Ask Jim about the cutoff calculation.

    /*
        int ncut = 1;// nspeicies**2

        int nspecies = parms->charmmParms->charmmPar->ljParmSize;

        int k = 0;
        for (int i = 0; i < nspecies; i++)
            for (int j = 0; j < nspecies; j++) {
                //double rcut = parms->rcut[j + nspecies * i].value;
                parms->rcut[k].value = parms->rmax;
                parms->rcut[k].mode = RCUT_ALL;
                parms->rcut[k].type = i * nspecies + j;
                k++;
            }
     */
    static RCUT_TYPE rcut[2];
    rcut[0].value = parms->rmax;
    rcut[0].mode = RCUT_LOCAL;
    rcut[0].type = -1;
    rcut[1].value = parms->rmax4all;
    rcut[1].mode = RCUT_ALL;
    rcut[1].type = -1;
    *n = 1;
    return rcut;
}

int getCHLJindexbyResAtmNameStd(char* resName, char* atmName, CHARMM_PARMS *charmmParms) //resClean
{
    LJCH_PARMS** ljParms = charmmParms->charmmPar->ljParms;
    int ljParmSize = charmmParms->charmmPar->ljParmSize;

    RESI_PARMS** resParms = charmmParms->charmmTop->resiParms;
    int resParmSize = charmmParms->charmmTop->resiParmSize;
    for (int i = 0; i < resParmSize; ++i)
    {
        if (resNameDiffer(resParms[i]->resName, resName) == 0)
        {
            TATOM_PARMS** atmList = resParms[i]->atomList;
            int atmListSize = resParms[i]->atomListSize;
            for (int j = 0; j < atmListSize; ++j)
            {
                if (strcmp(atmList[j]->atmName, atmName) == 0)
                {
                    for (int k = 0; k < ljParmSize; ++k)
                    {
                        if (strcmp(atmList[j]->atmTypePtr->atmType, ljParms[k]->atmTypeI) == 0) return k;
                    }
                }
            }
        }
    }

    return -1;
}

int getCHLJindexbyResAtmNameNon(char* resName, char*terName, char* atmName, CHARMM_PARMS *charmmParms) //resClean
{

    LJCH_PARMS** ljParms = charmmParms->charmmPar->ljParms;
    int ljParmSize = charmmParms->charmmPar->ljParmSize;

    RESI_PARMS** presParms = charmmParms->charmmTop->presParms;
    int presParmSize = charmmParms->charmmTop->presParmSize;
    for (int i = 0; i < presParmSize; ++i)
    {
        if (resNameDiffer(presParms[i]->resName, terName) == 0)
        {
            TATOM_PARMS** atmList = presParms[i]->atomList;
            int atmListSize = presParms[i]->atomListSize;
            for (int j = 0; j < atmListSize; ++j)
            {
                if (strcmp(atmList[j]->atmName, atmName) == 0)
                {
                    for (int k = 0; k < ljParmSize; ++k)
                    {
                        if (strcmp(atmList[j]->atmTypePtr->atmType, ljParms[k]->atmTypeI) == 0) return k;
                    }
                }
            }
        }
    }
    RESI_PARMS** resParms = charmmParms->charmmTop->resiParms;
    int resParmSize = charmmParms->charmmTop->resiParmSize;
    for (int i = 0; i < resParmSize; ++i)
    {
        if (resNameDiffer(resParms[i]->resName, resName) == 0)
        {
            TATOM_PARMS** atmList = resParms[i]->atomList;
            int atmListSize = resParms[i]->atomListSize;
            for (int j = 0; j < atmListSize; ++j)
            {
                if (strcmp(atmList[j]->atmName, atmName) == 0)
                {
                    for (int k = 0; k < ljParmSize; ++k)
                    {
                        if (strcmp(atmList[j]->atmTypePtr->atmType, ljParms[k]->atmTypeI) == 0) return k;
                    }
                }
            }
        }
    }

    return -1;
}

int getCHLJindexbySpecie(SPECIES* specie, CHARMM_PARMS *charmmParms)
{

    char *CTER = "CTER";
    char *NTER = "NTER";
    char *PROP = "PROP";
    char *GLYP = "GLYP";
    int index = 0;
    char resName[6]; //resClean

    char* specieStr = specie->name;

    char* delimPtr = strpbrk(specieStr, "cnx");
    char terminus = *delimPtr;
    int resNSize = delimPtr - specieStr;
    strncpy(resName, specieStr, resNSize); //resClean
    resName[resNSize] = '\0'; //resClean
    char *atmName = delimPtr + 1;

    if (terminus == 'x')
    {
        index = getCHLJindexbyResAtmNameStd(resName, atmName, charmmParms); //resClean
    }
    if (terminus == 'n')
    {
        char *terName = NTER;
        if (resNameDiffer(resName, "PRO") == 0) terName = PROP;
        if (resNameDiffer(resName, "GLYP") == 0) terName = GLYP;
        index = getCHLJindexbyResAtmNameNon(resName, terName, atmName, charmmParms); //resClean
    }
    if (terminus == 'c')
    {
        index = getCHLJindexbyResAtmNameNon(resName, CTER, atmName, charmmParms); //resClean
    }

    if (index < 0) printf("Wrong index for amino acid %s\n", specieStr);

    return index;
}

void charmmPair(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    profile(CHARMM_NONBOND, START);
    int IndexRmax = -1;
    for (int i = 0; i < sys->neighbor->nc; i++)
    {
        if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
    }
    if (IndexRmax == -1)
    {
        error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));
    }
    CHARMM_PARMS *charmmParms = parms->charmmParms;
    int nspecies = charmmParms->charmmPar->ljParmSize;
    NPARTICLE *particles = sys->neighbor->particles;
    double er = 0.0;
    double rcut = parms->rmax;

    int speciesIndex[sys->nspecies];
    for (int j = 0; j < sys->nspecies; j++)
        speciesIndex[sys->species[j]->index] = getCHLJindexbySpecie(sys->species[j], charmmParms);
    int sIndex[sys->nion];
    for (unsigned i = 0; i < sys->nion; i++)
        sIndex[i] = speciesIndex[sys->collection->state->species[i]->index];
    double* energy = sys->collection->state->potentialEnergy;
    //double* rx = sys->collection->state->rx;
    //double* ry = sys->collection->state->ry;
    //double* rz = sys->collection->state->rz;
    double* fx = sys->collection->state->fx;
    double* fy = sys->collection->state->fy;
    double* fz = sys->collection->state->fz;

    THREE_VECTOR fp;
    for (unsigned i = 0; i < sys->nlocal; i++)
    {
        int si = sIndex[i];
        PAIRS *pij = particles[i].ifirst[IndexRmax];
        while (pij != NULL)
        {
            RPAIR *p = pij->p;
            double r = p->r;
            int j = pij->j;
            int sj = sIndex[j];
            int sij = sj + nspecies*si;
            //double rcut = parms->rcut[sij].value;
            if (r < rcut)
            {
                //double dvdr;
                //double vij = parms->fcn(parms->parms[sij], r, &dvdr);
                // CHARMM LJ: V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
                CharmmLJ_PARMS *ljparms = parms->parms[sij];
                double sigma = ljparms->sigma;
                double eps = ljparms->eps;
                double ir = 1 / r;
                double sigma_r = sigma * ir;
                double s2 = sigma_r*sigma_r;
                double s4 = s2*s2;
                double s6 = s4*s2;
                double s12 = s6*s6;
                double dvdr = 12.0 * eps * (s6 - s12) * ir;
                double vij = eps * (s12 - 2.0 * s6);

                dvdr *= ir;

                er += vij;
                fp.x = -dvdr * p->x;
                fp.y = -dvdr * p->y;
                fp.z = -dvdr * p->z;

                fx[i] += fp.x;
                fy[i] += fp.y;
                fz[i] += fp.z;
                energy[i] += 0.5 * vij;

                fx[j] -= fp.x;
                fy[j] -= fp.y;
                fz[j] -= fp.z;
                energy[j] += 0.5 * vij;

            }
            pij = pij->ilink;
        }
    }

    /*
       double outer=units_convert(er, NULL, "kcal*mol^-1");
       printf("L-J Energy (kcal/mol) \n");
       printf("%15.5f \n", outer);

       printf("L-J Energy (default) \n");
       printf("%15.5f \n", er);     
     */
    parms->bioEnergies.nonBond += er;
    parms->bioEnergies.lj = er;
    e->eion += er;
    profile(CHARMM_NONBOND, END);
}

int compare(const void * a, const void * b)
{

    GID_ORDER *orderA = (GID_ORDER *) a;
    GID_ORDER *orderB = (GID_ORDER *) b;

    int d = 1; // use d a int type to avoid gid_type
    if (orderA->gid < orderB->gid)
    {
        d = -1;
    }

    return d;
}

void charmmGPU(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get current process id */
    profile(CHARMM_T, START);
    charmmConvalent(sys, parms, e); //Calculate all bonded terms
    profile(CHARMM_T, END);

    GPUCODE(charmmPairGPU(sys, parms, e);) //Subtract LJ energy within 1-4 atoms./

}



/*
   void pairqq0(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e) {
   double er, r, vij, dvdr, rcut;
   NPARTICLE *particles;
   PAIRS *pij;
   RPAIR *p;
   int IndexRmax = -1;
   for (int i = 0; i < sys->neighbor->nc; i++) {
   if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
   }
   if (IndexRmax == -1) {
   error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));
   }
   particles = sys->neighbor->particles;
   er = 0.0;
   rcut = parms->rcut[0].value;
   double *q = sys->collection->state->q;
   for (unsigned i = 0; i < sys->nlocal; i++) {
   int si = sys->collection->state->species[i]->index;
   pij = particles[i].ifirst[IndexRmax];
   while (pij != NULL) {
   p = pij->p;
   r = p->r;
   int j = pij->j;
   int sj = sys->collection->state->species[j]->index;
   int sij = sj + sys->nspecies*si;
   if (r < rcut) {
   double qij = q[i] * q[j];
   vij = qij * parms->fcn(parms->parms[sij], r, &dvdr);
   dvdr *= (qij / r);
   er += vij;
   p->e += vij;
   p->fp.x -= dvdr * p->x;
   p->fp.y -= dvdr * p->y;
   p->fp.z -= dvdr * p->z;
   }
   pij = pij->ilink;
   }
   }
   e->eion += er;
   }
 */

/*
 */

void charmm(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get current process id */
    profile(CHARMM_T, START);
    charmmConvalent(sys, parms, e); //Calculate all bonded terms
    charmmPair(sys, parms, e); //Subtract LJ energy within 1-4 atoms./
    printBioEnergies(parms->bioEnergies);
    profile(CHARMM_T, END);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */

