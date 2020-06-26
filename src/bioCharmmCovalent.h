#ifndef BIOCHARMMCOVALENT_H
#define BIOCHARMMCOVALENT_H
#include "state.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "energyInfo.h"
#include "three_algebra.h"

/* Single precision accuracy */
#define FLOAT_EPS    1e-08
#define NEAR_ZERO_ANGLE 0.017453292519943295 // 1 degree
#define NEAR_180_ANGLE  3.12413936106985     // 179 degree

typedef struct covalent_energy_struct
{
    double ebond;
    double eangle;
    double ecosangle;
    double eub;
    double etorsion;
    double eimpr;
    double ebpair;
    double eblj;
    double ebele;
    double ecmap;
} ECOVALENT;

#ifdef __cplusplus
extern "C"
{
#endif

void charmmResidues(SYSTEM*sys, CHARMMPOT_PARMS *parms);

#ifdef __cplusplus
}
#endif


int resNameDiffer(char* resName, char *name);
double resBondSorted(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resAngleSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resAngleCosineSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resAngleRestrainSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resUreyBradleySorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resTorsionSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resImproperSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resBondSortedWeighted(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resAngleSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resUreyBradleySortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resTorsionSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resImproperSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights);
double resCmap(STATE* state, CHARMM_PARMS* charmmParms, GID_ORDER* gidOrder, unsigned nlocal, unsigned ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e);
double resBpairSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, CHARMMPOT_PARMS *parms, ETYPE *e, double *eblj, double *ebele, BIOWEIGHTS* weights);
double resCGBpairSorted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, CHARMMPOT_PARMS *parms, ETYPE *e, double *eblj, double *ebele, BIOWEIGHTS* weights);

void bioVec(STATE* state, int atom1, int atom2, THREE_VECTOR* vec);
THREE_VECTOR bioScale(THREE_VECTOR vec, double scale);
double bioNorm(THREE_VECTOR vec);
double bioNorm2(THREE_VECTOR vec);
void bioUnit(STATE* state, int atom1, int atom2, THREE_VECTOR* unit);
double bioAngle(STATE* state, int atom1, int atom2, int atom3);
THREE_VECTOR bioVecAdd(double scale1, THREE_VECTOR vec1, double scale2, THREE_VECTOR vec2);
int bioDihedral(STATE* state, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                THREE_VECTOR* v_ij, THREE_VECTOR* v_kj, THREE_VECTOR* v_kl,
                THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL);
int bioDihedralFast(STATE* state, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                    THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL, THREE_SMATRIX* virial);

int resNameDiffer(char* resName, char *name);
#endif

