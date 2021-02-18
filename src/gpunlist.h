#ifndef GPUNLIST_H
#define GPUNLIST_H
#include "species.h"
#include "group.h"
#include "state.h"
#include "cudaUtils.h"

struct gputypes_st;

typedef struct gpuvirials_st {
    double *xx;
    double *yy;
    double *zz;
    double *xy;
    double *xz;
    double *yz;
} GPUVIRIALS;


typedef struct gpunlist_st {
    // binning data
    double *rxbg, *rybg, *rzbg; //gpu particle positions, permuted so particles are in bins
    double *results;
    double *e_all; //energy of each particle on gpu
    double *scratch, *scratch1, *scratch2, *scratch3, *scratch4; //scratch space
    double *partialCOMx, *partialCOMy, *partialCOMz;
    double *virCorx, *virCory, *virCorz;
    int *r_backg, *r_backbg, *binHeadsg, *binHeadsgSave, *binCountsg, *binCountsg2, *listIdsg, *partial_sumsg;
    double *minsg, *lensg;
    int *nbinsg, *nbrIds, *nbrIdsx, *nbrIdsy, *nbrIdsz;
    int numNbrs, nBinsTot;
    int *numNbrsxyz;
    unsigned *binsg, *binsbg;
    int allocFactor;
    int nMolecules;
    int *moleculeOffsets, *moleculeList;
    int *moleculeOffsetsg, *moleculeListg;
    int *molReferenceAtoms, *molReferenceAtomsg;
    double *molMasses, *molMassesg;

    gid_type * label;
    gid_type * labelSorted;
    // nlist data
    int **plist_h; //list of gpu pages, host array
    int **plist; //list of gpu pages, gpu array
    int *nbrSizes;
    int *mem;
    
    //move from the collection
    double *rx_gsort, *ry_gsort, *rz_gsort; //positions on gpu, sorted by gid
    double *rx_padg, *ry_padg, *rz_padg; //positions on gpu, sorted by gid, and "padded"

    double *e1; //energy of each particle on gpu
    double *e_all_h; //energy of each particle on gpu
    double *charge_g, *charge_bg; //gpu particle charges without bin permutation, and with bin permut.ation
    double *mass_h, *mass_g;
    int *species_h, *species_bg, *species_g;
    double *mmgpu;
    
    GPUVIRIALS *gpu_sion;
    GPUVIRIALS *gpu_sion_h;
    GPUVIRIALS *sion_h;

    int *pbc_g;
    THREE_MATRIX * hmat_g;
    THREE_MATRIX * hmati_g;

    //GPU random numbers
    double *gpuRandomx, *gpuRandomy, *gpuRandomz;
    //GPU special types
    struct gputypes_st *gpu_types;    
    
} GPUNLIST;

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
