#include <string.h>
#include <limits.h>
#include <math.h>
#include "system.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "energyInfo.h"
#include "expandbuffer.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "bioGid.h"
#include "preduce.h"
#include "codata.h"
#include "units.h"
#include "gid.h"
#include "three_algebra.h"


//#define RAD2DEG 57.29577951308232 // 180/PI
//#define BOH2ANG 0.52918 // Bohr to Angstrom

double bioBond(STATE* state, int atom1, int atom2)
{
    double bx = state->rx[atom1] - state->rx[atom2];
    double by = state->ry[atom1] - state->ry[atom2];
    double bz = state->rz[atom1] - state->rz[atom2];

    nearestImage(&bx, &by, &bz);
    //return sqrt(bx*bx+by*by+bz*bz)*BOH2ANG;
    return sqrt(bx * bx + by * by + bz * bz);
}

void bioVec(STATE* state, int atom1, int atom2, THREE_VECTOR* vec)
{
    //vec->x=-BOH2ANG*(state->rx[atom1]-state->rx[atom2]);
    //vec->y=-BOH2ANG*(state->ry[atom1]-state->ry[atom2]);
    //vec->z=-BOH2ANG*(state->rz[atom1]-state->rz[atom2]);    
    double bx = state->rx[atom1] - state->rx[atom2];
    double by = state->ry[atom1] - state->ry[atom2];
    double bz = state->rz[atom1] - state->rz[atom2];

    nearestImage(&bx, &by, &bz);

    vec->x = bx;
    vec->y = by;
    vec->z = bz;
}

double bioNorm(THREE_VECTOR vec)
{
    return (sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z));
}

double bioNorm2(THREE_VECTOR vec)
{
    return (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

THREE_VECTOR bioScale(THREE_VECTOR vec, double scale)
{
    vec.x *= scale;
    vec.y *= scale;
    vec.z *= scale;

    return vec;
}

THREE_VECTOR bioVecAdd(double scale1, THREE_VECTOR vec1, double scale2, THREE_VECTOR vec2)
{
    THREE_VECTOR vec;
    vec.x = scale1 * vec1.x + scale2 * vec2.x;
    vec.y = scale1 * vec1.y + scale2 * vec2.y;
    vec.z = scale1 * vec1.z + scale2 * vec2.z;
    return vec;
}

void bioUnit(STATE* state, int atom1, int atom2, THREE_VECTOR* unit)
{
    THREE_VECTOR vec;
    bioVec(state, atom1, atom2, &vec);
    double b = bioNorm(vec);

    unit->x = vec.x / b;
    unit->y = vec.y / b;
    unit->z = vec.z / b;
}

double bioAngle(STATE* state, int atom1, int atom2, int atom3)
{
    THREE_VECTOR unit21;
    bioUnit(state, atom2, atom1, &unit21);
    THREE_VECTOR unit23;
    bioUnit(state, atom2, atom3, &unit23);
    double angle = acos(unit21.x * unit23.x + unit21.y * unit23.y + unit21.z * unit23.z);

    return angle; // is rad not deg
}

double bioTorsion(STATE* state, int atom1, int atom2, int atom3, int atom4)
{
    THREE_VECTOR unit21;
    bioUnit(state, atom2, atom1, &unit21);
    THREE_VECTOR unit23;
    bioUnit(state, atom2, atom3, &unit23);

    double eabc_x = unit21.y * unit23.z - unit21.z * unit23.y;
    double eabc_y = unit21.z * unit23.x - unit21.x * unit23.z;
    double eabc_z = unit21.x * unit23.y - unit21.y * unit23.x;

    THREE_VECTOR unit32;
    bioUnit(state, atom3, atom2, &unit32);
    THREE_VECTOR unit34;
    bioUnit(state, atom3, atom4, &unit34);

    double ebcd_x = unit32.y * unit34.z - unit32.z * unit34.y;
    double ebcd_y = unit32.z * unit34.x - unit32.x * unit34.z;
    double ebcd_z = unit32.x * unit34.y - unit32.y * unit34.x;

    double exx = eabc_x * ebcd_x;
    double eyy = eabc_y * ebcd_y;
    double ezz = eabc_z * ebcd_z;

    double angle123 = bioAngle(state, atom1, atom2, atom3);
    double angle234 = bioAngle(state, atom2, atom3, atom4);
    double tau = (exx + eyy + ezz) / (sin(angle123) * sin(angle234));

    if (tau < -1.0) tau = acos(-1.0);
    else if (tau > 1.0) tau = acos(1.0);
    else tau = acos(tau);

    // Compute the sign of the torsion 
    double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
    double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
    double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
    double norm = cross_x * cross_x + cross_y * cross_y + cross_z*cross_z;
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;

    double sign = 1.0;
    double dot = cross_x * unit23.x + cross_y * unit23.y + cross_z * unit23.z;
    if (dot < 0.0) sign = -1.0;

    return tau*sign;
}

int bioDihedral(STATE* state, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                THREE_VECTOR* v_ij, THREE_VECTOR* v_kj, THREE_VECTOR* v_kl,
                THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL)
{
    //THREE_VECTOR v_ij;
    bioVec(state, indexI, indexJ, v_ij);
    //THREE_VECTOR v_kj;
    bioVec(state, indexK, indexJ, v_kj);
    //THREE_VECTOR v_kl;
    bioVec(state, indexK, indexL, v_kl);

    THREE_VECTOR v_A = cross(v_ij, v_kj);
    THREE_VECTOR v_B = cross(v_kj, v_kl);
    THREE_VECTOR v_C = cross(&v_A, v_kj);

    double v1_A = bioNorm(v_A);
    double v1_B = bioNorm(v_B);
    double v1_C = bioNorm(v_C);

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


    double dotCB = dot1(v_C, v_B);
    double sinX = dotCB / (v1_C * v1_B);
    *sinXptr = sinX;

    double dotAB = dot1(v_A, v_B);
    double cosX = dotAB / (v1_A * v1_B);
    //printf("cos=%f\n", cosX);

    /*
        THREE_VECTOR crossAB=cross(&v_A, &v_B);
        double dotkjAB=dot1(*v_kj, crossAB);
    
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
    (*dcosXdrI) = cross(&v_a, v_kj);

    // dcosX/drL = - r_kj X b   
    (*dcosXdrL) = cross(&v_b, v_kj);

    // dcosX/drJ = - r_ik X a + r_kl X b
    THREE_VECTOR v_ik;
    bioVec(state, indexI, indexK, &v_ik);
    THREE_VECTOR va_ik = cross(&v_ik, &v_a);
    THREE_VECTOR vb_kl = cross(v_kl, &v_b);
    (*dcosXdrJ) = bioVecAdd(-1, va_ik, 1, vb_kl);

    // dcosX/drK = - r_jl X b + r_ij X a
    THREE_VECTOR v_jl;
    bioVec(state, indexJ, indexL, &v_jl);
    THREE_VECTOR vb_jl = cross(&v_jl, &v_b);
    THREE_VECTOR va_ij = cross(v_ij, &v_a);
    (*dcosXdrK) = bioVecAdd(-1, vb_jl, 1, va_ij);

    return 0;


}

int bioDihedralFast(STATE* state, int indexI, int indexJ, int indexK, int indexL, double* angleX, double* sinXptr,
                    THREE_VECTOR* dcosXdrI, THREE_VECTOR* dcosXdrJ, THREE_VECTOR* dcosXdrK, THREE_VECTOR* dcosXdrL, THREE_SMATRIX* virial)
{

    double eps = 1e-12;
    THREE_VECTOR v_ij;
    THREE_VECTOR v_jk;
    THREE_VECTOR v_kl;
    bioVec(state, indexI, indexJ, &v_ij);
    bioVec(state, indexJ, indexK, &v_jk);
    bioVec(state, indexK, indexL, &v_kl);

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

    THREE_VECTOR m,n,mXn;
    CROSS(m, v_ij, v_jk);
    CROSS(n, v_jk, v_kl);
    CROSS(mXn, m, n); // Flip the conventional in Bekker 1995 JCC.
    double signnum = DOT(v_jk, mXn);
    double sign = (signnum < 0.0) ? -1.0 : 1.0;
    x = MAX(MIN(x, 1.0), -1.0);
    *angleX = sign * acos(x);
    *sinXptr = sin(*angleX);
/*
    CROSS(m, v_ij, v_jk);
    double dotX= b2*DOT(m,v_kl)-bc*DOT(m,v_jk);
    double signX = (dotX  < 0.0) ? -1.0 : 1.0;
    *sinXptr = signX*sqrt(1-x*x);
    *angleX = signX*acos(x); 
    printf("%e %e %e\n",dotX,*angleX,*sinXptr);
*/



    (*dcosXdrI) = ca;
    VOP2((*dcosXdrJ), =, cb, -, ca);
    VOP2((*dcosXdrK), =, cc, -, cb);
    VSVOP((*dcosXdrL), =, -1, *, cc);

    virial->xx = -1 * (ca.x * v_ij.x + cb.x * v_jk.x + cc.x * v_kl.x);
    virial->xy = -1 * (ca.x * v_ij.y + cb.x * v_jk.y + cc.x * v_kl.y);
    virial->xz = -1 * (ca.x * v_ij.z + cb.x * v_jk.z + cc.x * v_kl.z);
    virial->yy = -1 * (ca.y * v_ij.y + cb.y * v_jk.y + cc.y * v_kl.y);
    virial->yz = -1 * (ca.y * v_ij.z + cb.y * v_jk.z + cc.y * v_kl.z);
    virial->zz = -1 * (ca.z * v_ij.z + cb.z * v_jk.z + cc.z * v_kl.z);

   
   static FILE *file=NULL; 
   if (file == NULL) file = fopen("ddd.data","w"); 
    double cosab = ab/(sqrt(a2*b2));
    double cosbc = bc/(sqrt(c2*b2));
   
    return 0;

}


//Functions related to CMAP energy calculation

void determineCMAPCoef(double* x, double* coef)
{
    int size = 16;
    //store as 16by16 for easier management
    int Ainv[16][16] = {
        { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
        {-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0},
        { 9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1},
        {-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1},
        { 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
        {-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1},
        { 4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1}
    };

    //alpha = Ainv*xT 
    for (int i = 0; i < size; i++)
    {
        double temp = 0.0;
        for (int j = 0; j < size; j++)
        {
            temp += Ainv[i][j] * x[j];
        }
        coef[i] = temp;
    }
}

//return cmap energy
//Assumes that you have given it the correct CMAP PARM determined by cmapType
//Assumes you have calculated the relevant dihedrals phi, psi

double calculateCMAPValues(double phi, double psi, CMAP_PARMS* cmapParms, double *derivatives)
{
    int nGrid = 24;
    int size = 4;
    double resolution = 360.0 / ((double) nGrid);
    int res = (int) resolution;

    int indexPhi = (int) (phi / resolution);
    int indexPsi = (int) (psi / resolution);
    int phiPlus1 = (indexPhi + 1) % nGrid;
    int psiPlus1 = (indexPsi + 1) % nGrid;

    //build x for Coef v2
    double *x = ddcMalloc(size * size * sizeof (double));
    double *coef = ddcMalloc(size * size * sizeof (double));
    double *c = ddcMalloc(size * size * sizeof (double));
    //order x = {f:{(0,0) (1,0) (0,1) (1,1)} fx:same fy:same fxy}
    //f
    x[0] = cmapParms->grid[indexPhi * nGrid + indexPsi];
    x[1] = cmapParms->grid[phiPlus1 * nGrid + indexPsi];
    x[2] = cmapParms->grid[indexPhi * nGrid + psiPlus1];
    x[3] = cmapParms->grid[phiPlus1 * nGrid + psiPlus1];
    //fx= gridVal * resolution 
    x[4] = cmapParms->map_y1[indexPhi * nGrid + indexPsi] * res;
    x[5] = cmapParms->map_y1[phiPlus1 * nGrid + indexPsi] * res;
    x[6] = cmapParms->map_y1[indexPhi * nGrid + psiPlus1] * res;
    x[7] = cmapParms->map_y1[phiPlus1 * nGrid + psiPlus1] * res;
    //fy = gridVal*resolution
    x[8] = cmapParms->map_y2[indexPhi * nGrid + indexPsi] * res;
    x[9] = cmapParms->map_y2[phiPlus1 * nGrid + indexPsi] * res;
    x[10] = cmapParms->map_y2[indexPhi * nGrid + psiPlus1] * res;
    x[11] = cmapParms->map_y2[phiPlus1 * nGrid + psiPlus1] * res;
    //fxy = gridVal*resolution*resolution
    x[12] = cmapParms->map_y12[indexPhi * nGrid + indexPsi] * res*res;
    x[13] = cmapParms->map_y12[phiPlus1 * nGrid + indexPsi] * res*res;
    x[14] = cmapParms->map_y12[indexPhi * nGrid + psiPlus1] * res*res;
    x[15] = cmapParms->map_y12[phiPlus1 * nGrid + psiPlus1] * res*res;

    //Calculate Cijs using grid information
    determineCMAPCoef(x, coef);
    //reorder coef
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            c[count] = coef[j * size + i];
            count++;
        }
    }

    //Calculate Energy and Derivative Values
    double energy = 0.0;
    double dfx = 0.0;
    double dfy = 0.0;
    double ddfx = 0.0;
    double ddfy = 0.0;
    double ddfxy = 0.0;

    double temp1 = (phi - (double) (indexPhi) * resolution) / resolution;
    double temp2 = (psi - (double) (indexPsi) * resolution) / resolution;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double cval = c[i * 4 + j];
            energy += cval * pow(temp1, i) * pow(temp2, j);
            if (i != 0)
                dfx += cval * pow(temp1, i - 1) * pow(temp2, j) * i;
            if (j != 0)
                dfy += cval * pow(temp1, i) * pow(temp2, j - 1) * j;
            if (i != 0 && i != 1)
                ddfx += cval * pow(temp1, i - 2) * pow(temp2, j) * i * (i - 1);
            if (j != 0 && j != 1)
                ddfy += cval * pow(temp1, i) * pow(temp2, j - 2) * j * (j - 1);
            if (i != 0 && j != 0)
                ddfxy += cval * pow(temp1, i - 1) * pow(temp2, j - 1) * i * j;
        }
    }

    double deg2radian = 180.0 / M_PI;
    double rres = deg2radian / resolution;

    dfx *= rres;
    dfy *= rres;
    ddfx *= rres*rres;
    ddfy *= rres*rres;
    ddfxy *= rres*rres;

    derivatives[0] = units_convert(dfx, "kcal*mol^-1", NULL);
    derivatives[1] = units_convert(dfy, "kcal*mol^-1", NULL);
    derivatives[2] = units_convert(ddfx, "kcal*mol^-1", NULL);
    derivatives[3] = units_convert(ddfy, "kcal*mol^-1", NULL);
    derivatives[4] = units_convert(ddfxy, "kcal*mol^-1", NULL);

    //free
    ddcFree(x);
    ddcFree(coef);
    ddcFree(c);

    return units_convert(energy, "kcal*mol^-1", NULL);
}

void getPrevResSet(long long unsigned molID, long long unsigned resID, STATE* state, unsigned ntotal, gid_type* prevResSet, int *size)
{
    for (unsigned int i = 0; i < ntotal; i++)
    {
        gid_type currGid = state->label[i];
        long long unsigned currMolID = (currGid & molMask) >> 32;
        long long unsigned currResID = (currGid & resMask) >> 16;
        if (currMolID == molID && currResID == resID)
        {
            prevResSet[*size] = currGid;
            *size += 1;
        }
    }
}

gid_type getAtomGid(RESI_CONN* resiConn, char* atom, unsigned molID, unsigned resID)
{
    gid_type gid = -1;
    gid_type mID = molID;
    gid_type rID = resID;
    for (int i = 0; i < resiConn->groupListSize; i++)
    {
        TGROUP_PARMS* grpParms = resiConn->groupList[i];
        for (int j = 0; j < grpParms->grpAtomSize; j++)
        {
            TATOM_PARMS* atomParms = grpParms->grpAtoms[j];
            if (strncmp(atom, atomParms->atmName, 5) == 0)
            {
                gid_type groupID = grpParms->grpID;
                gid_type atomID = atomParms->atmID;
                gid = (mID << 32)+(rID << 16)+(groupID << 8)+(atomID << 8);
            }
        }
    }
    return gid;
}

double resCmap(STATE* state, CHARMM_PARMS* charmmParms, GID_ORDER* gidOrder, unsigned nlocal, unsigned ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e)
{

    double ectot = 0.0;
    if (resiConn->cmapListSize == 0) return ectot;
    CMAP_CONN **cmapList = resiConn->cmapList;
    for (int i = 0; i < resiConn->atomListSize; i++)
    {
        //for this atom, check if it is the N atom involved in the central cmap residue
        RANGE cmapRange = resiConn->atmRanges[i]->cmapRange;
        if (cmapRange.start == -1)
        {
            continue;
        }
        int cmapIndex = cmapRange.start;
        CMAP_CONN* cmapConn = cmapList[cmapIndex];
        int indexI = resRange->start + cmapConn->atmI;
        int indexJ = resRange->start + cmapConn->atmJ;
        int indexK = resRange->start + cmapConn->atmK;

        if ((unsigned) gidOrder[indexI].id >= nlocal)
        {
            continue;
        }
        //get central residue GID
        gid_type currGid = state->label[indexI];
        long long unsigned currMolID = (currGid & molMask) >> 32;
        long long unsigned currResID = (currGid & resMask) >> 16;
        long long unsigned nextResID = currResID + 1;
        gid_type nextGid = ((currMolID) << 32)+((nextResID) << 16)+((0) << 8)+((0) << 8);

        gid_type *val = (gid_type *) bsearch(&nextGid, state->label, ntotal, sizeof (gid_type), compareGid);

        if (val == NULL)
        {
            printf("CMAP: Error in Binary Search. Did not match gids!\n");
            exit(0);
        }
        unsigned nextIndex = ((gid_type *) val - (gid_type *) state->label);

        //determine CMAP type based on curr and next residue species 
        int prolineFlag = 0;
        char nextRes[5];
        strncpy(nextRes, state->species[nextIndex]->name, 4);

        char xDelim = 'x';
        char nDelim = 'n';
        char cDelim = 'c';

        char* speciesStr = state->species[nextIndex]->name;
        char* xDelimPtr = strchr(speciesStr, xDelim);
        char* nDelimPtr = strchr(speciesStr, nDelim);
        char* cDelimPtr = strchr(speciesStr, cDelim);

        char nextResName[5];

        if (xDelimPtr != NULL)
        {
            int resNSize = xDelimPtr - speciesStr;
            strncpy(nextResName, speciesStr, resNSize);
            nextResName[resNSize] = '\0';
        }
        else if (nDelimPtr != NULL)
        {
            int resNSize = nDelimPtr - speciesStr;
            strncpy(nextResName, speciesStr, resNSize);
            nextResName[resNSize] = '\0';
        }
        else if (cDelimPtr != NULL)
        {
            int resNSize = cDelimPtr - speciesStr;
            strncpy(nextResName, speciesStr, resNSize);
            nextResName[resNSize] = '\0';
        }

        if (strcmp(nextResName, "PRO") == 0)
        {
            prolineFlag = 1;
        }
        int cmapType = -1;
        cmapType = cmapList[prolineFlag]->cmapType;

        //previous res C info
        //This allocates more memory than needed 
        gid_type *prevRes = ddcMalloc(ntotal * sizeof (gid_type));
        int prevSize = 0;

        getPrevResSet(currMolID, currResID - 1, state, ntotal, prevRes, &prevSize);
        //Assumes both C and O are local 
        gid_type prevGid = -1;
        if (prevSize >= 2)
            prevGid = prevRes[prevSize - 2];
        else
            printf("Error in CMAP Term\n");

        gid_type *val2 = (gid_type *) bsearch(&prevGid, state->label, ntotal, sizeof (gid_type), compareGid);
        if (val2 == NULL)
        {
            printf("CMAP: Error in Binary Search. Did not match gids for previous Index!\n");
            exit(0);
        }
        unsigned prevIndex = ((gid_type *) val2 - (gid_type *) state->label);

        //calculate phi and psi
        double phi;
        double psi;
        double sinPhi;
        double sinPsi;
        //first dihedral
        //THREE_VECTOR v_ij;
        //THREE_VECTOR v_kj;
        //THREE_VECTOR v_kl;

        THREE_VECTOR dTdrI;
        THREE_VECTOR dTdrJ;
        THREE_VECTOR dTdrK;
        THREE_VECTOR dTdrL;

        THREE_SMATRIX virial1;

        bioDihedralFast(state, prevIndex, indexI, indexJ, indexK, &phi, &sinPhi, &dTdrI, &dTdrJ, &dTdrK, &dTdrL, &virial1);
        //second dihedral
        //THREE_VECTOR v2_ij;
        //THREE_VECTOR v2_kj;
        //THREE_VECTOR v2_kl;

        THREE_VECTOR dTdrI2;
        THREE_VECTOR dTdrJ2;
        THREE_VECTOR dTdrK2;
        THREE_VECTOR dTdrL2;

        THREE_SMATRIX virial2;

        bioDihedralFast(state, indexI, indexJ, indexK, nextIndex, &psi, &sinPsi, &dTdrI2, &dTdrJ2, &dTdrK2, &dTdrL2, &virial2);

        //Make sure compatible with CMAP grid definitions
        phi *= RAD2DEG;
        psi *= RAD2DEG;
        phi -= 180.0;
        psi -= 180.0;
        phi *= -1.0;
        psi *= -1.0;
        //use phi and psi for cmap calculation
        //TODO: Check if parallelization problem with this approach
        double * derivatives = ddcMalloc(5 * sizeof (double));
        CMAP_PARMS* cmapParms = charmmParms->charmmPar->cmapParms[cmapType];
        ectot = calculateCMAPValues(phi, psi, cmapParms, derivatives);
        //      for(int i=0; i<5; i++)
        //        printf("\tCHECK: derivatives[%d] = %f\n", i, derivatives[i]);

        // update energy / forces 

        //first dihedral
        THREE_VECTOR FdTdRi = bioScale(dTdrI, derivatives[0]);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, derivatives[0]);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, derivatives[0]);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, derivatives[0]);

        state->fx[prevIndex] -= FdTdRi.x;
        state->fy[prevIndex] -= FdTdRi.y;
        state->fz[prevIndex] -= FdTdRi.z;

        state->fx[indexI] -= FdTdRj.x;
        state->fy[indexI] -= FdTdRj.y;
        state->fz[indexI] -= FdTdRj.z;

        state->fx[indexJ] -= FdTdRk.x;
        state->fy[indexJ] -= FdTdRk.y;
        state->fz[indexJ] -= FdTdRk.z;

        state->fx[indexK] -= FdTdRl.x;
        state->fy[indexK] -= FdTdRl.y;
        state->fz[indexK] -= FdTdRl.z;

        // Pressure  

        e->virial.xx += virial1.xx * derivatives[0];
        e->virial.xy += virial1.xy * derivatives[0];
        e->virial.xz += virial1.xz * derivatives[0];
        e->virial.yy += virial1.yy * derivatives[0];
        e->virial.yz += virial1.yz * derivatives[0];
        e->virial.zz += virial1.zz * derivatives[0];

        //second dihedral
        THREE_VECTOR FdTdRi2 = bioScale(dTdrI2, derivatives[1]);
        THREE_VECTOR FdTdRj2 = bioScale(dTdrJ2, derivatives[1]);
        THREE_VECTOR FdTdRk2 = bioScale(dTdrK2, derivatives[1]);
        THREE_VECTOR FdTdRl2 = bioScale(dTdrL2, derivatives[1]);

        state->fx[indexI] -= FdTdRi2.x;
        state->fy[indexI] -= FdTdRi2.y;
        state->fz[indexI] -= FdTdRi2.z;

        state->fx[indexJ] -= FdTdRj2.x;
        state->fy[indexJ] -= FdTdRj2.y;
        state->fz[indexJ] -= FdTdRj2.z;

        state->fx[indexK] -= FdTdRk2.x;
        state->fy[indexK] -= FdTdRk2.y;
        state->fz[indexK] -= FdTdRk2.z;

        state->fx[nextIndex] -= FdTdRl2.x;
        state->fy[nextIndex] -= FdTdRl2.y;
        state->fz[nextIndex] -= FdTdRl2.z;

        // Pressure

        e->virial.xx += virial2.xx * derivatives[1];
        e->virial.xy += virial2.xy * derivatives[1];
        e->virial.xz += virial2.xz * derivatives[1];
        e->virial.yy += virial2.yy * derivatives[1];
        e->virial.yz += virial2.yz * derivatives[1];
        e->virial.zz += virial2.zz * derivatives[1];

    }
    return ectot;
}

void connectiveEnergy(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms, ETYPE *e)
{

    CHARMM_PARMS* charmmParms = parms->charmmParms;
    RESI_CONN* resiConn = findResiConnNew(charmmParms, name);


    if (parms->bond_fcn)
        parms->bioEnergies.bond += parms->bond_fcn(state, gidOrder, nlocal, nTotal, resRange, resiConn, e, charmmParms->charmmWeights); //BOND 

    if (parms->angle_fcn)
        parms->bioEnergies.angle += parms->angle_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); //Angle

    if (parms->cosangle_fcn)
        parms->bioEnergies.angle += parms->cosangle_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); // Cosine-based Angle

    if (parms->rebangle_fcn)
        parms->bioEnergies.angle += parms->rebangle_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); // Cosine-based Angle
    
    if (parms->ureybradley_fcn)
        parms->bioEnergies.ub += parms->ureybradley_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); //Urey-Bradley

    if (parms->torsion_fcn)
        parms->bioEnergies.torsion += parms->torsion_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); //Torsion

    if (parms->improper_fcn)
        parms->bioEnergies.impr += parms->improper_fcn(state, gidOrder, nlocal, resRange, resiConn, e, charmmParms->charmmWeights); //Improper

    if (parms->cmap_fcn)
        parms->bioEnergies.cmap += parms->cmap_fcn(state, charmmParms, gidOrder, nlocal, nTotal, resRange, resiConn, e); //CMAP

    if (parms->bpair_fcn)
    {
        double eblj = 0;
        double ebele = 0;
        double bpair = parms->bpair_fcn(state, gidOrder, nlocal, resRange, resiConn, parms, e, &eblj, &ebele, charmmParms->charmmWeights); //BondPair
        parms->bioEnergies.bpair += bpair;
        parms->bioEnergies.blj += eblj;
        parms->bioEnergies.bele += ebele;
    }

}


/* Local Variables: */
/* tab-width: 3 */
/* End: */

