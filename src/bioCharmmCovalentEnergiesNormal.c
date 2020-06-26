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

double resBondNormal(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double ebtot = 0; //
    for (int i = 0; i < resiConn->bondListSize; ++i)
    {
        BOND_CONN* bondConn = resiConn->bondList[i];
        int indexI = resRange->start + bondConn->atmI;
        int indexJ = resRange->start + bondConn->atmJ;
        // Energy
        //double b=bioBond(state, indexI, indexJ);
        THREE_VECTOR bVec;
        bioVec(state, indexI, indexJ, &bVec);
        double b = bioNorm(bVec);
        //printf("name %s  %s\n", state->species[indexI]->name,state->species[indexJ]->name);
        BOND_PARMS* bondPtr = bondConn->bondPtr;
        double bDelta = b - bondPtr->b0;
        double eb = bondPtr->kb * bDelta*bDelta;
        ebtot = ebtot + eb;
        if (eb > 1.0)
        {
            //FIX ME PUT THESE 3 PRINTS BACK
            /*
            printf("I=%d J=%d SpecieI=%s SpecieJ=%s GidI=%llu, Gid=%llu, BondI=%s BondJ=%s\n", 
                    indexI, indexJ, state->species[indexI]->name, state->species[indexJ]->name,
                    state->label[indexI], state->label[indexJ],
                    bondPtr->atmTypeI, bondPtr->atmTypeJ);         
            printf("Bond length=%f bondPara=%f bond energy=%f\n",b, bondPtr->b0, eb);
            printf("x1=%f    y1=%f    z1=%f    x2=%f    y2=%f    z2=%f\n", 
                    state->rx[indexI],state->ry[indexI],state->rz[indexI],
                    state->rx[indexJ],state->ry[indexJ],state->rz[indexJ]);
             */
        }
        // Force
        THREE_VECTOR unit_ij;
        //bioUnit(state, indexI, indexJ, &unit_ij);
        unit_ij.x = bVec.x / b;
        unit_ij.y = bVec.y / b;
        unit_ij.z = bVec.z / b;

        double fxDelta = 2 * bondPtr->kb * bDelta * unit_ij.x;
        double fyDelta = 2 * bondPtr->kb * bDelta * unit_ij.y;
        double fzDelta = 2 * bondPtr->kb * bDelta * unit_ij.z;
        state->fx[indexI] += fxDelta; // F=-dU/dx
        state->fy[indexI] += fyDelta;
        state->fz[indexI] += fzDelta;
        state->fx[indexJ] -= fxDelta; // Fj=-Fi
        state->fy[indexJ] -= fyDelta;
        state->fz[indexJ] -= fzDelta;

        // Pressure

        double fpxx = fxDelta * bVec.x;
        double fpxy = fxDelta * bVec.y;
        double fpxz = fxDelta * bVec.z;
        double fpyy = fyDelta * bVec.y;
        double fpyz = fyDelta * bVec.z;
        double fpzz = fzDelta * bVec.z;

        e->virial.xx += fpxx;
        e->virial.xy += fpxy;
        e->virial.xz += fpxz;
        e->virial.yy += fpyy;
        e->virial.yz += fpyz;
        e->virial.zz += fpzz;

        /*
                printf("Bond Delta force: x=%f y=%f z=%f\n",fxDelta, fyDelta, fzDelta);
                printf("Bond atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
                printf("Bond atomJ force: x=%f y=%f z=%f\n",state->fx[indexJ], state->fy[indexJ], state->fz[indexJ]);
         */
    }
    return ebtot;
}

double resAngleNormal(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eatot = 0;
    for (int i = 0; i < resiConn->angleListSize; ++i)
    {
        ANGLE_CONN* angleConn = resiConn->angleList[i];
        int indexI = resRange->start + angleConn->atmI;
        int indexJ = resRange->start + angleConn->atmJ;
        int indexK = resRange->start + angleConn->atmK;
        //printf("angle %s %s %s \n", state->species[indexI]->name,state->species[indexJ]->name, state->species[indexK]->name);
        // Energy
        //double a=bioAngle(state, indexI, indexJ, indexK);
        THREE_VECTOR vec_ij;
        bioVec(state, indexI, indexJ, &vec_ij);
        double b_ij = bioNorm(vec_ij);
        THREE_VECTOR unit_ij;
        unit_ij.x = vec_ij.x / b_ij;
        unit_ij.y = vec_ij.y / b_ij;
        unit_ij.z = vec_ij.z / b_ij;

        THREE_VECTOR vec_kj;
        bioVec(state, indexK, indexJ, &vec_kj);
        double b_kj = bioNorm(vec_kj);
        THREE_VECTOR unit_kj;
        unit_kj.x = vec_kj.x / b_kj;
        unit_kj.y = vec_kj.y / b_kj;
        unit_kj.z = vec_kj.z / b_kj;

        double a = acos(unit_ij.x * unit_kj.x + unit_ij.y * unit_kj.y + unit_ij.z * unit_kj.z);

        ANGLE_PARMS* anglePtr = angleConn->anglePtr;
        //double aDelta=a-anglePtr->theta0/RAD2DEG;
        double aDelta = a - anglePtr->theta0;
        double ea = anglePtr->ktheta * aDelta*aDelta; // convert deg to rad
        eatot = eatot + ea;
        if (ea > 0.5)
        {
            printf("I=%d J=%d K=%d\n", indexI, indexJ, indexK);
            printf("Angle =%f Angle Deg=%f AnglePara=%f Angle energy=%f\n", a, RAD2DEG*a, anglePtr->theta0, ea);
        }
        // Force
        /*
               eprintf("Angle atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
                printf("Angle atomJ force: x=%f y=%f z=%f\n",state->fx[indexJ], state->fy[indexJ], state->fz[indexJ]);
                printf("Angle atomK force: x=%f y=%f z=%f\n",state->fx[indexK], state->fy[indexK], state->fz[indexK]);
         */

        double sinabs = sin(a);
        //if(sinabs<0) sinabs=-sinabs;
        double coef_i = -2 * anglePtr->ktheta * aDelta / (b_ij * sinabs);
        double cosTheta = cos(a);
        double fxDeltaI = coef_i * (unit_kj.x - unit_ij.x * cosTheta);
        double fyDeltaI = coef_i * (unit_kj.y - unit_ij.y * cosTheta);
        double fzDeltaI = coef_i * (unit_kj.z - unit_ij.z * cosTheta);

        double coef_k = -2 * anglePtr->ktheta * aDelta / (b_kj * sinabs);
        double fxDeltaK = coef_k * (unit_ij.x - unit_kj.x * cosTheta);
        double fyDeltaK = coef_k * (unit_ij.y - unit_kj.y * cosTheta);
        double fzDeltaK = coef_k * (unit_ij.z - unit_kj.z * cosTheta);

        state->fx[indexI] += fxDeltaI;
        state->fy[indexI] += fyDeltaI;
        state->fz[indexI] += fzDeltaI;

        state->fx[indexK] += fxDeltaK;
        state->fy[indexK] += fyDeltaK;
        state->fz[indexK] += fzDeltaK;

        state->fx[indexJ] -= (fxDeltaI + fxDeltaK);
        state->fy[indexJ] -= (fyDeltaI + fyDeltaK);
        state->fz[indexJ] -= (fzDeltaI + fzDeltaK);

        // Pressure

        double fpxx = (fxDeltaI * vec_ij.x + fxDeltaK * vec_kj.x);
        double fpxy = (fxDeltaI * vec_ij.y + fxDeltaK * vec_kj.y);
        double fpxz = (fxDeltaI * vec_ij.z + fxDeltaK * vec_kj.z);
        double fpyy = (fyDeltaI * vec_ij.y + fyDeltaK * vec_kj.y);
        double fpyz = (fyDeltaI * vec_ij.z + fyDeltaK * vec_kj.z);
        double fpzz = (fzDeltaI * vec_ij.z + fzDeltaK * vec_kj.z);

        e->virial.xx += fpxx;
        e->virial.xy += fpxy;
        e->virial.xz += fpxz;
        e->virial.yy += fpyy;
        e->virial.yz += fpyz;
        e->virial.zz += fpzz;

        /*
                printf("Angle Delta force I: x=%f y=%f z=%f\n",fxDeltaI, fyDeltaI, fzDeltaI);
                printf("Angle Delta force K: x=%f y=%f z=%f\n",fxDeltaK, fyDeltaK, fzDeltaK);
        
                printf("Angle atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
                printf("Angle atomJ force: x=%f y=%f z=%f\n",state->fx[indexJ], state->fy[indexJ], state->fz[indexJ]);
                printf("Angle atomK force: x=%f y=%f z=%f\n",state->fx[indexK], state->fy[indexK], state->fz[indexK]);
         */

    }
    return eatot;
}

double resUreyBradleyNormal(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eubtot = 0;
    for (int i = 0; i < resiConn->angleListSize; ++i)
    {
        ANGLE_CONN* angleConn = resiConn->angleList[i];
        int indexI = resRange->start + angleConn->atmI;
        int indexK = resRange->start + angleConn->atmK;
        ANGLE_PARMS* anglePtr = angleConn->anglePtr;
        if (anglePtr->kub > 0)
        {
            // Energy
            //double bub=bioBond(state, indexI, indexK);            
            THREE_VECTOR bVec;
            bioVec(state, indexI, indexK, &bVec);
            double bub = bioNorm(bVec);

            double bubDelta = bub - anglePtr->s0;
            double eub = anglePtr->kub * bubDelta*bubDelta;
            eubtot = eubtot + eub;
            if (eub > 0.1)
            {
                printf("I=%d K=%d \n", indexI, indexK);
                printf("Urey-Bradley 1-3 bond=%f Ub s0=%f  UB energy=%f\n", bub, anglePtr->s0, eub);
            }
            // Force
            THREE_VECTOR unit_ik;
            //bioUnit(state, indexI, indexK, &unit_ik);
            unit_ik.x = bVec.x / bub;
            unit_ik.y = bVec.y / bub;
            unit_ik.z = bVec.z / bub;

            double fxDelta = 2 * anglePtr->kub * bubDelta * unit_ik.x;
            double fyDelta = 2 * anglePtr->kub * bubDelta * unit_ik.y;
            double fzDelta = 2 * anglePtr->kub * bubDelta * unit_ik.z;
            state->fx[indexI] += fxDelta;
            state->fy[indexI] += fyDelta;
            state->fz[indexI] += fzDelta;
            state->fx[indexK] -= fxDelta;
            state->fy[indexK] -= fyDelta;
            state->fz[indexK] -= fzDelta;

            // Pressure

            double fpxx = fxDelta * bVec.x;
            double fpxy = fxDelta * bVec.y;
            double fpxz = fxDelta * bVec.z;
            double fpyy = fyDelta * bVec.y;
            double fpyz = fyDelta * bVec.z;
            double fpzz = fzDelta * bVec.z;

            e->virial.xx += fpxx;
            e->virial.xy += fpxy;
            e->virial.xz += fpxz;
            e->virial.yy += fpyy;
            e->virial.yz += fpyz;
            e->virial.zz += fpzz;

            /*
                        printf("Urey-Bradley Delta force: x=%f y=%f z=%f\n",fxDelta, fyDelta, fzDelta);
                        printf("Urey-Bradley atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
                        printf("Urey-Bradley atomJ force: x=%f y=%f z=%f\n",state->fx[indexK], state->fy[indexK], state->fz[indexK]);            
             */
        }
    }
    return eubtot;
}

double resTorsionNormal(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double etorstot = 0;
    for (int i = 0; i < resiConn->torsListSize; ++i)
    {
        TORS_CONN* torsConn = resiConn->torsList[i];
        int indexI = resRange->start + torsConn->atmI;
        int indexJ = resRange->start + torsConn->atmJ;
        int indexK = resRange->start + torsConn->atmK;
        int indexL = resRange->start + torsConn->atmL;
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

        if (bioDihedral(state, indexI, indexJ, indexK, indexL, &tors, &sinX,
                        &v_ij, &v_kj, &v_kl, &dTdrI, &dTdrJ, &dTdrK, &dTdrL) < 0)
        {
            continue; // skip the torsion with angle approach to 0.
            // also skip the force calculation because a & b approach to 0;

        }

        TORSION_PARMS* torsPtr = torsConn->torsPtr;
        double etors = torsPtr->kchi * (1 + cos(torsPtr->n * tors - torsPtr->delta));
        etorstot = etorstot + etors;

        //double tors2=bioTorsion(state, indexI, indexJ, indexK, indexL);
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
            ktors = torsPtr->kchi * torsPtr->n * sin(torsPtr->n * tors - torsPtr->delta) / sinX;
            //negative sign canceled out
        }
        else
        {
            double nX = torsPtr->n*tors;
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
            double ratio = torsPtr->n * (1 - nX2 / 6 + nX4 / 120 - nX6 / 5040 + nX8 / 362880 - nX10 / 39916800) / (1 - X2 / 6 + X4 / 120 - X6 / 5040 + X8 / 362880 - X10 / 39916800);
            if (torsPtr->delta < NEAR_ZERO_ANGLE)
            { // for delta=0
                ktors = torsPtr->kchi * torsPtr->n*ratio;
            }
            else if (torsPtr->delta > NEAR_180_ANGLE)
            { // for delta=180
                ktors = -torsPtr->kchi * torsPtr->n*ratio;
                // sin(x-pi)=-sin(pi-x)=-sin(x)
            }
            else
            {
                //printf("Warning: torsPtr->delta doesn't equal to 0 or 180");
                ktors = torsPtr->kchi * torsPtr->n*ratio;
            }

        }
        //double ktors=-torsPtr->kchi*torsPtr->n*sin(torsPtr->n*tors-torsPtr->delta);

        THREE_VECTOR FdTdRi = bioScale(dTdrI, ktors);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, ktors);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, ktors);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, ktors);

        state->fx[indexI] -= FdTdRi.x;
        state->fy[indexI] -= FdTdRi.y;
        state->fz[indexI] -= FdTdRi.z;

        state->fx[indexJ] -= FdTdRj.x;
        state->fy[indexJ] -= FdTdRj.y;
        state->fz[indexJ] -= FdTdRj.z;

        state->fx[indexK] -= FdTdRk.x;
        state->fy[indexK] -= FdTdRk.y;
        state->fz[indexK] -= FdTdRk.z;

        state->fx[indexL] -= FdTdRl.x;
        state->fy[indexL] -= FdTdRl.y;
        state->fz[indexL] -= FdTdRl.z;

        // Pressure

        double fpxx = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.x + (FdTdRk.x - FdTdRj.x) * v_kj.x + (FdTdRk.x - FdTdRl.x) * v_kl.x); // v_kl=-v_lk
        double fpxy = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.y + (FdTdRk.x - FdTdRj.x) * v_kj.y + (FdTdRk.x - FdTdRl.x) * v_kl.y);
        double fpxz = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.z + (FdTdRk.x - FdTdRj.x) * v_kj.z + (FdTdRk.x - FdTdRl.x) * v_kl.z);
        double fpyy = 0.5 * ((FdTdRi.y - FdTdRj.y) * v_ij.y + (FdTdRk.y - FdTdRj.y) * v_kj.y + (FdTdRk.y - FdTdRl.y) * v_kl.y);
        double fpyz = 0.5 * ((FdTdRi.y - FdTdRj.y) * v_ij.z + (FdTdRk.y - FdTdRj.y) * v_kj.z + (FdTdRk.y - FdTdRl.y) * v_kl.z);
        double fpzz = 0.5 * ((FdTdRi.z - FdTdRj.z) * v_ij.z + (FdTdRk.z - FdTdRj.z) * v_kj.z + (FdTdRk.z - FdTdRl.z) * v_kl.z);

        e->virial.xx += fpxx;
        e->virial.xy += fpxy;
        e->virial.xz += fpxz;
        e->virial.yy += fpyy;
        e->virial.yz += fpyz;
        e->virial.zz += fpzz;
    }

    return etorstot;
}

double resImproperNormal(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eimprtot = 0;
    for (int i = 0; i < resiConn->imprListSize; ++i)
    {
        IMPR_CONN* imprConn = resiConn->imprList[i];
        int indexI = resRange->start + imprConn->atmI; // indexI is the center atom
        int indexJ = resRange->start + imprConn->atmJ;
        int indexK = resRange->start + imprConn->atmK;
        int indexL = resRange->start + imprConn->atmL;
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

        if (bioDihedral(state, indexI, indexJ, indexK, indexL, &impr, &sinX,
                        &v_ij, &v_kj, &v_kl, &dTdrI, &dTdrJ, &dTdrK, &dTdrL) < 0)
        {
            continue; // skip the torsion with angle approach to 0.
            // also skip the force calculation because a & b approach to 0;
        }

        IMPROPER_PARMS* imprPtr = imprConn->imprPtr;
        //double imprDelta=impr-imprPtr->psi0/RAD2DEG;
        double imprDelta = impr - imprPtr->psi0;
        double eimpr = imprPtr->kpsi * imprDelta*imprDelta;
        eimprtot = eimprtot + eimpr;
        if (eimpr > 0.1)
        {
            printf("I=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            printf("Improper =%f Improper angle=%f kpsi=%f psi0=%f Improper Energy=%f\n",
                   impr, impr*RAD2DEG, imprPtr->kpsi, imprPtr->psi0*RAD2DEG, eimpr);
        }
        //Force
        //Improper is the permutation of Torsion.

        double absX = sinX;
        if (sinX < 0) absX = -sinX;

        double kimpr;

        if (absX > FLOAT_EPS)
        {
            kimpr = -2 * imprPtr->kpsi * imprDelta / sinX;
        }
        else
        {
            if (imprPtr->psi0 > NEAR_ZERO_ANGLE)
            {
                printf("Warning: imprPtr->kpsi doesn't equal to 0");
            }
            double impr2 = impr*impr;
            double impr4 = impr2*impr2;
            double impr6 = impr4*impr2;
            double impr8 = impr4*impr4;
            double impr10 = impr8*impr2;
            kimpr = -2 * imprPtr->kpsi / (1 - impr2 / 6 + impr4 / 120 - impr6 / 5040 + impr8 / 362880 - impr10 / 39916800);
        }

        //double kimpr=2*imprPtr->kpsi*imprDelta;

        THREE_VECTOR FdTdRi = bioScale(dTdrI, kimpr);
        THREE_VECTOR FdTdRj = bioScale(dTdrJ, kimpr);
        THREE_VECTOR FdTdRk = bioScale(dTdrK, kimpr);
        THREE_VECTOR FdTdRl = bioScale(dTdrL, kimpr);

        state->fx[indexI] -= FdTdRi.x;
        state->fy[indexI] -= FdTdRi.y;
        state->fz[indexI] -= FdTdRi.z;

        state->fx[indexJ] -= FdTdRj.x;
        state->fy[indexJ] -= FdTdRj.y;
        state->fz[indexJ] -= FdTdRj.z;

        state->fx[indexK] -= FdTdRk.x;
        state->fy[indexK] -= FdTdRk.y;
        state->fz[indexK] -= FdTdRk.z;

        state->fx[indexL] -= FdTdRl.x;
        state->fy[indexL] -= FdTdRl.y;
        state->fz[indexL] -= FdTdRl.z;

        // Pressure

        double fpxx = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.x + (FdTdRk.x - FdTdRj.x) * v_kj.x + (FdTdRk.x - FdTdRl.x) * v_kl.x); // v_kl=-v_lk
        double fpxy = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.y + (FdTdRk.x - FdTdRj.x) * v_kj.y + (FdTdRk.x - FdTdRl.x) * v_kl.y);
        double fpxz = 0.5 * ((FdTdRi.x - FdTdRj.x) * v_ij.z + (FdTdRk.x - FdTdRj.x) * v_kj.z + (FdTdRk.x - FdTdRl.x) * v_kl.z);
        double fpyy = 0.5 * ((FdTdRi.y - FdTdRj.y) * v_ij.y + (FdTdRk.y - FdTdRj.y) * v_kj.y + (FdTdRk.y - FdTdRl.y) * v_kl.y);
        double fpyz = 0.5 * ((FdTdRi.y - FdTdRj.y) * v_ij.z + (FdTdRk.y - FdTdRj.y) * v_kj.z + (FdTdRk.y - FdTdRl.y) * v_kl.z);
        double fpzz = 0.5 * ((FdTdRi.z - FdTdRj.z) * v_ij.z + (FdTdRk.z - FdTdRj.z) * v_kj.z + (FdTdRk.z - FdTdRl.z) * v_kl.z);

        e->virial.xx += fpxx;
        e->virial.xy += fpxy;
        e->virial.xz += fpxz;
        e->virial.yy += fpyy;
        e->virial.yz += fpyz;
        e->virial.zz += fpzz;

    }

    return eimprtot;
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */

