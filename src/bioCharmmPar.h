#ifndef BIOCHARMMPAR_H
#define BIOCHARMMPAR_H

//BONDS
//V(bond) = Kb(b - b0)**2
//Kb: kcal/mole/A**2
//b0: Angstrom

//atom type Kb          b0
//NH2   CT1   240.00      1.455  ! From LSN NH2-CT2

typedef struct bond_parms_str
{
    char atmTypeI[5];
    char atmTypeJ[5];
    double kb;
    double b0;
} BOND_PARMS;

//ANGLES
//V(angle) = Ktheta(Theta - Theta0)**2
//V(Urey-Bradley) = Kub(S - S0)**2
//Ktheta: kcal/mole/rad**2
//Theta0: degrees
//Kub: kcal/mole/A**2 (Urey-Bradley)
//S0: Angstrom

//atom types     Ktheta    Theta0   Kub     S0
//H    NH2  CT1   50.000    111.00              ! From LSN HC-NH2-CT2
//NH2  CT1  HB    38.000    109.50   50.00   2.1400 ! From LSN NH2-CT2-HA

typedef struct angle_parms_str
{
    char atmTypeI[5];
    char atmTypeJ[5];
    char atmTypeK[5];
    double ktheta;
    double theta0;
    double kub;
    double s0;
} ANGLE_PARMS;


//DIHEDRALS
//V(dihedral) = Kchi(1 + cos(n(chi) - delta))
//Kchi: kcal/mole
//n: multiplicity
//delta: degrees

//atom types             Kchi    n   delta
//NH2  CT1  C    O        0.0000  1     0.00

typedef struct torsion_parms_str
{
    char atmTypeI[5];
    char atmTypeJ[5];
    char atmTypeK[5];
    char atmTypeL[5];
    int n;
    double kchi;
    double delta;
} TORSION_PARMS;

//IMPROPER
//V(improper) = Kpsi(psi - psi0)**2
//Kpsi: kcal/mole/rad**2
//psi0: degrees
//note that the second column of numbers (0) is ignored

//atom types           Kpsi                   psi0
//CPB  CPA  NPH  CPA    20.8000         0      0.0000 ! ALLOW HEM

typedef struct improper_parms_str
{
    char atmTypeI[5];
    char atmTypeJ[5];
    char atmTypeK[5];
    char atmTypeL[5];
    double kpsi;
    double psi0;
} IMPROPER_PARMS;

//CMAP

//NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
//cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
//!adm jr., 5/08/91, suggested cutoff scheme

//V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
//epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
//Rmin/2: Angstrom, Rmin,i,j = Rmin/2,i + Rmin/2,j

//atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
//C      0.000000  -0.110000     2.000000 ! ALLOW   PEP POL ARO
//CP1    0.000000  -0.020000     2.275000   0.000000  -0.010000     1.900000 ! ALLOW   AL 

typedef struct ljch_parms_str
{
    char atmTypeI[5];
    double epsilon;
    double rmin;
    double eps14;
    double rmin14;
} LJCH_PARMS;

//SARA: add a point to CMAP struct here
//It is still unclear to me how data is stored for Charmm Parms
//This struct is likely not consistent with the other structs
//Assumption: We have 6 different CMAP types. We will create
//an array(?, list?) to hold each type. This struct describes
//the general information needed for each CMAP type.

typedef struct cmap_parms_str
{
    double grid[24 * 24];
    double map_y1[24 * 24];
    double map_y2[24 * 24];
    double map_y12[24 * 24];
} CMAP_PARMS;

typedef struct charmm_par_str
{
    int bondParmSize;
    int angleParmSize;
    int torParmSize;
    int imprParmsize;
    int ljParmSize;
    int cmapParmSize;

    BOND_PARMS ** bondParms;
    ANGLE_PARMS ** angleParms;
    TORSION_PARMS **torParms;
    IMPROPER_PARMS **imprParms;
    LJCH_PARMS **ljParms;
    CMAP_PARMS **cmapParms;
} CHARMM_PAR;

int parseCharmmPar(const char *fileName, CHARMM_PAR *charmmPars);

int freeCharmmPar(CHARMM_PAR *charmmPars);

#endif /* BIOCHARMMPAR_H */

