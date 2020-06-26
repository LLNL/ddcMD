#ifndef EWALD_H
#define EWALD_H
#include <complex.h>
#include "system.h"
#include "energyInfo.h"
#include "three_algebra.h"
#include "threadScheduler.h"

#define ewaldV(u2)    (M_2_SQRTPI*(-1.0 + (u2)/3.0 - (u2)*(u2)/10.0))
#define ewaldDVDrOverR(u2)  (M_2_SQRTPI*( 2.0/3.0 - (u2)*4.0/10.0)); 
enum KSPACE_CLASS { NO_KSPACE, EWALD_KSPACE, PPPM_KSPACE };
enum RSPACE_CLASS { NO_RSPACE, EWALD_RSPACE };
enum EWALD_TERMS  { RSPACE_TERM=1, KSPACE_TERM=2, MADELUNG_TERM=4,INFINITY_TERM=8};
enum CHARGE_TYPE { CONSTANT_CHARGE, VARIABLE_CHARGE }; 
typedef struct qdipole_st { double q,q2,q4;THREE_VECTOR dipole;} QDIPOLE ; 
typedef struct rspace_parms_str 
{ 
    double alpha, kappa, q, q2, q4, rho, vol;
    double  rspaceError, rmax;  
    double psiBar; 
    double psiBarStar; 
    double gamma; 
    FOUR_VECTOR chi; 
    THREE_SMATRIX  virialF, virialE; 
    double (*rspaceFcn)(void *,double,double *);
    THREAD_SCHEDULER *schedule; 
} RSPACE_PARMS; 

typedef struct ewald_parms_str
{
    double alpha, kappa, q, q2, q4, rho, vol; 
    FOUR_VECTOR chi; 
    THREE_MATRIX h; 
    double error,rspaceError,kspaceError;
    THREE_VECTOR dipole; 
    double epsilon; 
    double psi; 
    double er,ek,eM,ed;
    THREE_SMATRIX  virialR,virialK,virialM,virialD;
    int chargeCalcFlag; 
    int chargeType; 
    double cpu_r, cpu_k;
    int optimize_alpha; 
    int rspace_type ; 
    int kspace_type ; 
    RSPACE_PARMS *rspace_parms; 
    void *kspace_parms; 
    FOUR_VECTOR *farField; 
    FOUR_VECTOR *nearField; 
    int fieldWriteRate; 
    char *filename; 
    int ewaldPotentialIndex; 
    int evaluate; 
} EWALD_PARMS;

typedef double ( *EWALD_RFUNCTION) ( RSPACE_PARMS *parms , double r, double *dvdr_over_r); 

double ewaldErr(double u, double k, double *dErr);
void ewaldBegin(SYSTEM *sys,EWALD_PARMS *parms,ETYPE *energyInfo) ;
void ewaldEnd(SYSTEM *sys,EWALD_PARMS *parms,ETYPE *energyInfo) ;
double ewald_get_rmax(EWALD_PARMS *parms);
void ewald_put_rspaceError(EWALD_PARMS *parms, double rspaceError);
double ewald_rspaceError(SYSTEM *sys, RSPACE_PARMS *parms,double r);
double chargeFromSpecies(SYSTEM *sys,double *q1, double *q2) ;

typedef struct kspace_parms_str 
{ 
    double alpha, kappa, q, q2, q4, rho, vol;
    FOUR_VECTOR chi ; 
    //THREE_MATRIX  h; 
    THREE_SMATRIX virial; 
    double  kspaceError, kmax; 
    double wnSum; 
    int nkvectors, nmax; 
    THREE_INT *nv; 
    double complex *px,*py,*pz;
    double complex  *vqppp,*vqppm,*vqpmp,*vqpmm;
} EWALD_KSPACE_PARMS; 

void ewald_kspace_parms( 
        EWALD_KSPACE_PARMS *parms, 
        double kspaceError, 
        double alpha, 
        double kappa, 
        double q, 
        double q2, 
        double q4, 
        double rho, 
        THREE_MATRIX h);


void printEwaldInfo(FILE *file, SYSTEM *sys, EWALD_PARMS *parms);
double ewald_kspaceError(SYSTEM *sys, void  *kspace_parms,double kmax);
EWALD_KSPACE_PARMS *ewald_kspaceInit();
void ewald_kspaceFree(EWALD_KSPACE_PARMS *parms);
double psiBar(double u, double xi);
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
