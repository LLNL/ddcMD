#ifndef BIOCHARMMPARMS_H
#define BIOCHARMMPARMS_H

#include "potential.h"
#include "state.h"
//#include "bioCharmm.h"
#include "bioCharmmPar.h"
#include "bioCharmmTop.h"
#include "bioTransform.h"
#include "gid.h"

#define HYDROGENMASSCUTOFF 2.35 //Use 1.1 for hydrogen mass cutoff

#define ATOMMAXBONDNUM 32
#define INTERMASSPNUM 4

/*
typedef struct charmmlj_struct 
{ 
    double rcut,eps,sigma;
    double shift;
}  CharmmLJ_PARMS;
 */

enum EXCLUDEPOTENIALTERMS
{
    bondMask = 1, angleMask = 2, cosangleMask = 4, ureybradleyMask = 8, torsionMask = 16, improperMask = 32, cMapMask = 64, nonBondMask = 128, rebangleMask = 256
};

typedef struct specInfo_struct
{
    int rIndex, gIndex, aIndex, offset, isResidueOwner;
    char *name;
} SPECIESINFO;

typedef struct charmmlj_gpu_struct
{
    double *rcut, *eps, *sigma;
    double *shift;
    double *charge;
    int *sIndex;
    int *cg_species_index;
    int *cg_species_index_b;
    int nspecies;
    double ke;
    double crf;
    double krf;
    double iepsilon_r;
    /* //
    int *angleListI, *angleListJ, *angleListK;
    int *torsListI, *torsListJ, *torsListK, *torsListL;
     */
} CharmmLJGPU_PARMS;

typedef struct charmmbond_gpu_struct
{
    double *rxs, *rys, *rzs;
    gid_type * labelSorted;

    THREE_MATRIX *h1;
    int * r_backBondInp;
    int * r_backBond;

    int * resID;
    int * resIDs;

    int *numAtomsInResidue;
    int *numAtomsInResiduePrefixSum;
    int atomListSize;
    //bond list parms
    int bondListSize;
    int *resiBondListStarts, *resiBondListSizes;
    int *atmBondListStarts, *atmBondListSizes;
    int *bondListI, *bondListJ;
    double *kb;
    double *b0;
} CharmmBONDGPU_PARMS;

typedef struct constraint_gpu_struct
{
    int numPair, constraintListSize, numAtom;
    int *resiConstraintCount; //size=#types(constraints)
    int *resiConstraintCountPrefixSum;
    int *constraintAtoms;
    int *constraintAtomsPrefixSum;
    int *constraintAtomIDs;
    int * constraintPairs; //size =#types(constraint pairs)
    int * constraintPairsPrefixSum;
    int * constraintPairsI;
    int * constraintPairsJ;
    int * constraintPairsIindex;
    int * constraintPairsJindex;
    double *constraintPairsDistance;
} CONSTRAINTGPU_PARMS;

typedef struct bpair_gpu_struct
{
    int bpairListSize;
    int *resiBpairListStarts, *resiBpairListSizes;
    int *atmBpairListStarts, *atmBpairListSizes;
    int *bpairListI, *bpairListJ;
    double *sigma;
    double *eps;
    double *shift;
    double *q;
} CharmmBPAIRGPU_PARMS;

typedef struct charmmangle_gpu_struct
{
    double *ktheta;
    double *theta0;
    double *kub;
    double *s0;

    //angle list parms
    int angleListSize;
    int *resiAngleListStarts, *resiAngleListSizes;
    int *atmAngleListStarts, *atmAngleListSizes;
    int *angleListI, *angleListJ, *angleListK;

} CharmmANGLEGPU_PARMS;

typedef struct charmmtorsion_gpu_struct
{
    double *kchi;
    double *delta;
    int *n;
    //torsion list parms
    int torsionListSize;
    int *resiTorsionListStarts, *resiTorsionListSizes;
    int *atmTorsionListStarts, *atmTorsionListSizes;
    int *torsionListI, *torsionListJ, *torsionListK, *torsionListL;

} CharmmTORSIONGPU_PARMS;

typedef struct charmmimproper_gpu_struct
{
    double *kpsi;
    double *psi0;

    //improper list parms
    int improperListSize;
    int *resiImproperListStarts, *resiImproperListSizes;
    int *atmImproperListStarts, *atmImproperListSizes;
    int *improperListI, *improperListJ, *improperListK, *improperListL;

} CharmmIMPROPERGPU_PARMS;

typedef struct bond_conn_str
{
    int atmI;
    int atmJ;
    BOND_PARMS *bondPtr;
} BOND_CONN;

typedef struct angle_conn_str
{
    int atmI;
    int atmJ;
    int atmK;
    ANGLE_PARMS *anglePtr;
} ANGLE_CONN;

typedef struct bndcnt_list_str
{
    int atom[ATOMMAXBONDNUM]; // Assume atom has less than 6 bonds
} BNDCNT_LIST;

typedef struct torsion_conn_str
{
    int atmI;
    int atmJ;
    int atmK;
    int atmL;

    //    char atmTypeI[5];
    //    char atmTypeL[5];
    TORSION_PARMS *torsPtr;
} TORS_CONN;

typedef struct impr_conn_str
{
    int atmI;
    int atmJ;
    int atmK;
    int atmL;

    IMPROPER_PARMS *imprPtr;
} IMPR_CONN;

typedef struct bpair_conn_str
{
    int atmI;
    int atmJ;
    double sigma;
    double eps;
    double shift;

} BPAIR_CONN;

typedef struct cmap_conn_str
{
    //central residue only
    int atmI;
    int atmJ;
    int atmK;
    //what type of grid parameters
    int cmapType;

    CMAP_PARMS *cmapPtr;
} CMAP_CONN;

typedef struct weight_conn_str
{
    //This struct is for atom
    //pairs to exclude from
    //charge calculations
    int atmI;
    int atmJ;
} ZEROWEIGHTS;

typedef struct range_str
{
    int start;
    int end;
} RANGE;

typedef struct atmrange_str
{
    RANGE bondRange;
    RANGE angleRange;
    RANGE cosangleRange;
    RANGE rebangleRange;
    RANGE torsRange;
    RANGE imprRange;
    RANGE bpairRange;
    RANGE cmapRange;
    RANGE weightRange;

} ATMRANGE;

typedef struct heavyatm_str
{
    int atomID;
    int group;
    int atom;

    int numHydrogen;
    // Assume less than 5 hydrogens bonded to heavy atom
    int hydrogen[4];
    int hGroup[4];
    int hGrpAtm[4];
    BOND_CONN* hBndConn[4];
} HEAVYATM;

typedef struct cons_pair_str
{
    int atomI;
    int atomIindex;
    int atomJ;
    int atomJindex;
    double distance;
} CONS_PAIR;

typedef struct constraint_str
{
    int atomID;
    int numPair;
    int numAtom;
    int* atomIDList;
    CONS_PAIR** conspairList;

} CONSTRAINT;

typedef struct resi_conn_str
{
    int resID;
    int resType;
    int nTer;
    int cTer;
    char resName[5];
    double charge;

    int centerAtom;

    int atomListSize;
    int groupListSize;
    int bondListSize;
    int angleListSize;
    int cosangleListSize;
    int rebangleListSize;
    int torsListSize;
    int imprListSize;
    int bpairListSize;
    int cmapListSize;
    int weightListSize;

    int heavyAtmsSize;
    int consListSize;

    MASS_PARMS* nh1MassParm;
    MASS_PARMS* cMassParm;

    TATOM_PARMS ** atomList;
    TGROUP_PARMS ** groupList;

    int *firstBond;
    int *bondCnt;
    BNDCNT_LIST *bndCntList;
    BOND_CONN** bondList;
    ANGLE_CONN** angleList;
    ANGLE_CONN** cosangleList;
    ANGLE_CONN** rebangleList;
    TORS_CONN** torsList;
    IMPR_CONN** imprList;
    BPAIR_CONN** bpairList;
    CMAP_CONN** cmapList;

    ATMRANGE** atmRanges;

    HEAVYATM** heavyAtms;
    CONSTRAINT** consList;

    ZEROWEIGHTS** weightList;
} RESI_CONN;

typedef struct list_node_struct
{
    gid_type resiStart;
    unsigned resID;
    int size;
    int first;
    int cnt;
    char * name;
    RESI_CONN* resiConn;
} LISTNODE;

typedef struct set_list_struct
{
    int molSize;
    int listSize;
    LISTNODE *list;
} SETLIST;

typedef struct
{
    gid_type gid;
    int id;
} GID_ORDER;

typedef struct charmm_parms_str
{
    POTENTIAL *potential;
    unsigned nspecies;
    void **parms;

    CHARMM_TOP* charmmTop;
    CHARMM_PAR* charmmPar;

    int resiConnSize;
    RESI_CONN** resiConnList;

    //For Molecule Transformations
    BIOWEIGHTS* charmmWeights;

} CHARMM_PARMS;

typedef struct bioEnergy_struct
{
    double bond;
    double angle;
    double ub;
    double torsion;
    double impr;
    double cmap;
    double blj;
    double bele;
    double bpair;
    double lj;
    double ele;
    double nonBond;
} BIOENERGIES;

typedef struct charmpot_parms_str
{
    POTENTIAL *potential;
    unsigned nspecies;
    char **species_list;
    //   RCUT_TYPE *rcut;
    double rmax;
    double rmax4all;
    // For martini reaction field parameters
    double rcoulomb;
    double epsilon_r;
    double epsilon_rf;
    double krf;
    double crf;
    // End of martini reaction field parameters
    BIOENERGIES bioEnergies;

    double(*nonBond_fcn)(void *, double, double *);
    void **parms;

    //device structs
    CharmmLJGPU_PARMS* gpu_ljparms;
    CharmmBONDGPU_PARMS *gpu_bondparms;
    CharmmANGLEGPU_PARMS *gpu_angleparms;
    CharmmANGLEGPU_PARMS *gpu_cosangleparms;
    CharmmANGLEGPU_PARMS *gpu_rebangleparms;
    CharmmTORSIONGPU_PARMS *gpu_torsionparms;
    CharmmIMPROPERGPU_PARMS *gpu_improperparms;
    CharmmBPAIRGPU_PARMS *gpu_bpairparms;
    CONSTRAINTGPU_PARMS * gpu_constraintparms;
    //host structs for temporary holding of serialized data
    CharmmBONDGPU_PARMS *bondparms_h;
    CharmmANGLEGPU_PARMS *angleparms_h;
    CharmmANGLEGPU_PARMS *cosangleparms_h;
    CharmmANGLEGPU_PARMS *rebangleparms_h;
    CharmmTORSIONGPU_PARMS *torsionparms_h;
    CharmmIMPROPERGPU_PARMS *improperparms_h;
    CharmmBPAIRGPU_PARMS *bpairparms_h;
    CONSTRAINTGPU_PARMS * constraintparms_h;
    //host structs with pointers to gpu mem
    CharmmLJGPU_PARMS* gpu_ljparms_h;
    CharmmBONDGPU_PARMS *gpu_bondparms_h;
    CharmmANGLEGPU_PARMS *gpu_angleparms_h;
    CharmmANGLEGPU_PARMS *gpu_cosangleparms_h;
    CharmmANGLEGPU_PARMS *gpu_rebangleparms_h;
    CharmmTORSIONGPU_PARMS *gpu_torsionparms_h;
    CharmmIMPROPERGPU_PARMS *gpu_improperparms_h;
    CharmmBPAIRGPU_PARMS *gpu_bpairparms_h;
    CONSTRAINTGPU_PARMS *gpu_constraintparms_h;

    CHARMM_PARMS* charmmParms;
    int useAutoResidueOwnership; //determines if (a task) owning species 0 in a residue -> task owns entire residue.  
    int * isResidueOwnerSpecies; //this array is not used if above flag is set to 0 (false  
    SPECIESINFO *speciesInfo;
    GID_ORDER *gidOrder;
    SETLIST residueSet;
    STATE statechpad;
    int weightedCharmm;
    // Function pointers for connectivity energy calculation
    double(*bond_fcn)(void *, void *, int, int, void *, void *, void *, void *);
    double(*angle_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*cosangle_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*rebangle_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*ureybradley_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*torsion_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*improper_fcn)(void *, void *, int, void *, void *, void *, void *);
    double(*bpair_fcn)(void *, void *, int, void*, void *, void *, void *, double*, double*, void *);
    double(*cmap_fcn)(void *, void *, void *, int, int, void *, void *, void *);
    int excludePotentialTerm;

} CHARMMPOT_PARMS;

//static int numOfBonds;
//static int numOfAngles;
//static int numOfDihedral;
//static int numOfimporper;

int parseCharmmParms(const char *topFile, const char *parFile, CHARMM_PARMS* charmmParms);

#ifdef __cplusplus
extern "C"
{
#endif
RESI_CONN* findResiConn(CHARMM_PARMS* charmmParms, char* resName, int nTer, int cTer);
RESI_CONN* findResiConnNew(CHARMM_PARMS* charmmParms, char* resName);

RESI_CONN* findResiConnBySpecieName(CHARMM_PARMS* charmmParms, char* spName);
#ifdef __cplusplus
}
#endif
void sortBondList(CHARMM_PARMS*charmmParms);
void genAtmRange(CHARMM_PARMS* charmmParms);

#endif /* BIOCHARMMPARMS_H */

