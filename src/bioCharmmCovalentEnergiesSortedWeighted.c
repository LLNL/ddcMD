#include <string.h>
#include <limits.h>
#include <math.h>
#include "system.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "bioTransform.h"
#include "energyInfo.h"
#include "expandbuffer.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "bioGid.h"
#include "preduce.h"
#include "codata.h"
#include "units.h"

double resBondSortedWeighted(STATE* state, GID_ORDER* gidOrder, int nlocal, int ntotal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double ebtot = 0;
    BOND_CONN ** bondList = resiConn->bondList;
    for (int i = 0; i < resiConn->atomListSize; i++)
    {
        //for this atom, get the range of adjacent bonds where this atom is listed first in pair data structure
        RANGE bondRange = resiConn->atmRanges[i]->bondRange;
        if (bondRange.start == -1)
        {
            continue;
        } // -1 means no relevant bonds...so go to next pair
        //bondRanges are the intra-residue offsets of each atom in the bond pair...use to get the global indeces 
        for (int bondIndex = bondRange.start; bondIndex < bondRange.end; bondIndex++)
        {
            BOND_CONN* bondConn = bondList[bondIndex];
            int indexI = resRange->start + (bondConn->atmI);
            int indexJ = resRange->start + (bondConn->atmJ);
            if (indexI >= ntotal || indexJ >= ntotal)
            {
                continue;
            } //skip if atom's don't exist
            if (gidOrder[indexI].id >= (int) nlocal)
            {
                continue;
            }

            // numOfBonds++;

            // Energy
            //double b=bioBond(state, indexI, indexJ);
            THREE_VECTOR bVec;
            bioVec(state, indexI, indexJ, &bVec);
            double b = bioNorm(bVec);

            BOND_PARMS* bondPtr = bondConn->bondPtr;
            double bDelta = b - bondPtr->b0;
            double eb = bondPtr->kb * bDelta*bDelta;
            //The problem is that the atomList is only the total number of atoms
            //types in the topology file, not the atomList of the system
            //double w = get2Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ]);
            double w = get2Weights(weights, resiConn->atomList[bondConn->atmI], resiConn->atomList[bondConn->atmJ]);
            eb *= w;
            ebtot = ebtot + eb;
            if (eb > 1.0)
            {

                printf("I=%d J=%d SpecieI=%s SpecieJ=%s GidI=%lu, Gid=%lu, BondI=%s BondJ=%s\n",
                       indexI, indexJ, state->species[indexI]->name, state->species[indexJ]->name,
                       state->label[indexI], state->label[indexJ],
                       bondPtr->atmTypeI, bondPtr->atmTypeJ);
                printf("Bond length=%f bondParaL=%f bondParaP=%f bond energy=%f\n",
                       units_convert(b, NULL, "Angstrom"),
                       units_convert(bondPtr->b0, NULL, "Angstrom"),
                       units_convert(bondPtr->kb, NULL, "kcal*mol^-1*Angstrom^-2"),
                       units_convert(eb, NULL, "kcal*mol^-1"));
                printf("x1=%f    y1=%f    z1=%f    x2=%f    y2=%f    z2=%f\n",
                       state->rx[indexI], state->ry[indexI], state->rz[indexI],
                       state->rx[indexJ], state->ry[indexJ], state->rz[indexJ]);

            }
            // Force
            THREE_VECTOR unit_ij;
            //bioUnit(state, indexI, indexJ, &unit_ij);
            unit_ij.x = bVec.x / b;
            unit_ij.y = bVec.y / b;
            unit_ij.z = bVec.z / b;

            double kforce = -2 * bondPtr->kb*bDelta;
            double fxDelta = w * kforce * unit_ij.x;
            double fyDelta = w * kforce * unit_ij.y;
            double fzDelta = w * kforce * unit_ij.z;
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
    }
    return ebtot;
}

double resAngleSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eatot = 0;
    ANGLE_CONN ** angleList = resiConn->angleList;
    //DEBUG ONLY
    int angleCount = 0;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE angleRange = resiConn->atmRanges[i]->angleRange;
        if (angleRange.start == -1)
        {
            continue;
        }
        for (int angleIndex = angleRange.start; angleIndex < angleRange.end; angleIndex++)
        {
            ANGLE_CONN* angleConn = angleList[angleIndex];
            int indexI = resRange->start + angleConn->atmI;
            int indexJ = resRange->start + angleConn->atmJ;
            int indexK = resRange->start + angleConn->atmK;

            //if not local atom than exit out
            if (gidOrder[indexJ].id >= (int) nlocal)
            {
                continue;
            }
            // Energy
            //double a=bioAngle(state, indexI, indexJ, indexK);
            //numOfAngles++;

            //printf("Angle: %s - %s - %s\n", state->species[indexI]->name, state->species[indexJ]->name, state->species[indexK]->name);

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

            double cosTheta = unit_ij.x * unit_kj.x + unit_ij.y * unit_kj.y + unit_ij.z * unit_kj.z;
            double a = acos(cosTheta);

            ANGLE_PARMS* anglePtr = angleConn->anglePtr;
            //double aDelta=a-anglePtr->theta0/RAD2DEG;
            double aDelta = a - anglePtr->theta0;
            double ea = anglePtr->ktheta * aDelta*aDelta; // convert deg to rad
            //double w = get3Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ], resiConn->atomList[indexK]);

            double w = get3Weights(weights, resiConn->atomList[angleConn->atmI], resiConn->atomList[angleConn->atmJ], resiConn->atomList[angleConn->atmK]);
            ea *= w;
            eatot = eatot + ea;

            angleCount++;

            if (ea > 0.5)
            {
                printf("I=%d J=%d K=%d\n", indexI, indexJ, indexK);
                printf("Angle =%f Angle Deg=%f AnglePara=%f Angle energy=%f\n", a, RAD2DEG*a, anglePtr->theta0, ea);
            }
            // Force
            /*
                    printf("Angle atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
                    printf("Angle atomJ force: x=%f y=%f z=%f\n",state->fx[indexJ], state->fy[indexJ], state->fz[indexJ]);
                    printf("Angle atomK force: x=%f y=%f z=%f\n",state->fx[indexK], state->fy[indexK], state->fz[indexK]);
             */

            double sinabs = sin(a);
            double coef_i = 2 * anglePtr->ktheta * aDelta / (b_ij * sinabs);

            double fxDeltaI = w * coef_i * (unit_kj.x - unit_ij.x * cosTheta);
            double fyDeltaI = w * coef_i * (unit_kj.y - unit_ij.y * cosTheta);
            double fzDeltaI = w * coef_i * (unit_kj.z - unit_ij.z * cosTheta);

            double coef_k = 2 * anglePtr->ktheta * aDelta / (b_kj * sinabs);
            double fxDeltaK = w * coef_k * (unit_ij.x - unit_kj.x * cosTheta);
            double fyDeltaK = w * coef_k * (unit_ij.y - unit_kj.y * cosTheta);
            double fzDeltaK = w * coef_k * (unit_ij.z - unit_kj.z * cosTheta);

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

    }
    return eatot;
}

double resAngleCosineSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eatot = 0;
    ANGLE_CONN ** cosangleList = resiConn->cosangleList;
    //DEBUG ONLY
    int angleCount = 0;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE cosangleRange = resiConn->atmRanges[i]->cosangleRange;
        if (cosangleRange.start == -1)
        {
            continue;
        }
        for (int angleIndex = cosangleRange.start; angleIndex < cosangleRange.end; angleIndex++)
        {
            ANGLE_CONN* angleConn = cosangleList[angleIndex];
            int indexI = resRange->start + angleConn->atmI;
            int indexJ = resRange->start + angleConn->atmJ;
            int indexK = resRange->start + angleConn->atmK;

            //if not local atom than exit out
            if (gidOrder[indexJ].id >= (int) nlocal)
            {
                continue;
            }
            // Energy
            //double a=bioAngle(state, indexI, indexJ, indexK);
            //numOfAngles++;

            //printf("Angle: %s - %s - %s\n", state->species[indexI]->name, state->species[indexJ]->name, state->species[indexK]->name);

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

            double cosA = unit_ij.x * unit_kj.x + unit_ij.y * unit_kj.y + unit_ij.z * unit_kj.z;

            ANGLE_PARMS* anglePtr = angleConn->anglePtr;
            double aDelta = cosA - anglePtr->theta0;
            double w = get3Weights(weights, resiConn->atomList[angleConn->atmI], resiConn->atomList[angleConn->atmJ], resiConn->atomList[angleConn->atmK]);
            double ea = anglePtr->ktheta * aDelta*aDelta; // convert deg to rad
            ea *= w;
            eatot = eatot + ea;

            angleCount++;

            if (ea > 0.5)
            {
                printf("I=%d J=%d K=%d\n", indexI, indexJ, indexK);
                printf("CosA =%f  AnglePara=%f Angle energy=%f\n", cosA, anglePtr->theta0, ea);

                printf("CosA =%f  AnglePara=%f Angle energy=%f\n", cosA, anglePtr->theta0, ea);
            }
            // Force
            /*
             printf("Angle atomI force: x=%f y=%f z=%f\n",state->fx[indexI], state->fy[indexI], state->fz[indexI]);
             printf("Angle atomJ force: x=%f y=%f z=%f\n",state->fx[indexJ], state->fy[indexJ], state->fz[indexJ]);
             printf("Angle atomK force: x=%f y=%f z=%f\n",state->fx[indexK], state->fy[indexK], state->fz[indexK]);
             */

            double coef_i = -2 * anglePtr->ktheta * aDelta / b_ij;
            double fxDeltaI = w * coef_i * (unit_kj.x - unit_ij.x * cosA);
            double fyDeltaI = w * coef_i * (unit_kj.y - unit_ij.y * cosA);
            double fzDeltaI = w * coef_i * (unit_kj.z - unit_ij.z * cosA);

            double coef_k = -2 * anglePtr->ktheta * aDelta / b_kj;
            double fxDeltaK = w * coef_k * (unit_ij.x - unit_kj.x * cosA);
            double fyDeltaK = w * coef_k * (unit_ij.y - unit_kj.y * cosA);
            double fzDeltaK = w * coef_k * (unit_ij.z - unit_kj.z * cosA);

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

    }
    //printf("CHECK: Angle count = %d\n", angleCount);
    return eatot;
}

double resUreyBradleySortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eubtot = 0;
    ANGLE_CONN ** angleList = resiConn->angleList;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE angleRange = resiConn->atmRanges[i]->angleRange;
        if (angleRange.start == -1)
        {
            continue;
        }
        for (int angleIndex = angleRange.start; angleIndex < angleRange.end; angleIndex++)
        {
            ANGLE_CONN* angleConn = angleList[angleIndex];

            int indexI = resRange->start + angleConn->atmI;
            int indexJ = resRange->start + angleConn->atmJ;
            int indexK = resRange->start + angleConn->atmK;

            //if not local atom than exit out
            if (gidOrder[indexJ].id >= (int) nlocal)
            {
                continue;
            }

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
                // double w = get3Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ], resiConn->atomList[indexK]);

                double w = get3Weights(weights, resiConn->atomList[angleConn->atmI], resiConn->atomList[angleConn->atmJ], resiConn->atomList[angleConn->atmK]);
                eub *= w;
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

                double kforce = -2 * anglePtr->kub*bubDelta;
                double fxDelta = w * kforce * unit_ik.x;
                double fyDelta = w * kforce * unit_ik.y;
                double fzDelta = w * kforce * unit_ik.z;
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
    }
    return eubtot;
}

double resTorsionSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{

    double etorstot = 0;
    TORS_CONN ** torsList = resiConn->torsList;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE torsRange = resiConn->atmRanges[i]->torsRange;
        if (torsRange.start == -1)
        {
            continue;
        }
        for (int torsIndex = torsRange.start; torsIndex < torsRange.end; torsIndex++)
        {
            TORS_CONN* torsConn = torsList[torsIndex];
            int indexI = resRange->start + torsConn->atmI;
            int indexJ = resRange->start + torsConn->atmJ;
            int indexK = resRange->start + torsConn->atmK;
            int indexL = resRange->start + torsConn->atmL;

            //if (( state->rx[indexI]==-1 || state->rx[indexJ]==-1 || state->rx[indexK]==-1 || state->rx[indexL]==-1)){continue;}   
            //if not local atom than exit out
            if (gidOrder[indexJ].id >= (int) nlocal)
            {
                continue;
            }

            //numOfDihedral++;
            //printf("Dihedral: %s - %s - %s - %s\n", 
            //        state->species[indexI]->name, state->species[indexJ]->name, state->species[indexK]->name, state->species[indexL]->name);

            //Energy
            double tors;
            double sinX;
            //THREE_VECTOR v_ij;
            //THREE_VECTOR v_kj;
            //THREE_VECTOR v_kl;

            THREE_VECTOR dTdrI;
            THREE_VECTOR dTdrJ;
            THREE_VECTOR dTdrK;
            THREE_VECTOR dTdrL;

            THREE_SMATRIX virial;

            if (bioDihedralFast(state, indexI, indexJ, indexK, indexL, &tors, &sinX,
                                &dTdrI, &dTdrJ, &dTdrK, &dTdrL, &virial) < 0)
            {
                continue; // skip the torsion with angle approach to 0.
                // also skip the force calculation because a & b approach to 0;

            }

            TORSION_PARMS* torsPtr = torsConn->torsPtr;
            double etors = torsPtr->kchi * (1 + cos(torsPtr->n * tors - torsPtr->delta));
            // double w = get4Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ], resiConn->atomList[indexK], resiConn->atomList[indexL]);

            double w = get4Weights(weights, resiConn->atomList[torsConn->atmI], resiConn->atomList[torsConn->atmJ], resiConn->atomList[torsConn->atmK], resiConn->atomList[torsConn->atmL]);
            etors *= w;
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
            //check dihedral
            // printf("CHECK Dihedral: residue = %s ID = %d nTer = %d cTer = %d\n", resiConn->resName, resiConn->resID, resiConn->cTer, resiConn->nTer); //resClean
            //  printf("\tI=%d J=%d K=%d L=%d\n", indexI, indexJ, indexK, indexL);
            //  printf("\tktors=%f eTors=%f sinx = %f\n", ktors, etors, sinX);
            //ktors=-ktors; // Flip the sign for bioVec

            THREE_VECTOR tempI = bioScale(dTdrI, ktors);
            THREE_VECTOR tempJ = bioScale(dTdrJ, ktors);
            THREE_VECTOR tempK = bioScale(dTdrK, ktors);
            THREE_VECTOR tempL = bioScale(dTdrL, ktors);

            THREE_VECTOR FdTdRi = bioScale(tempI, w);
            THREE_VECTOR FdTdRj = bioScale(tempJ, w);
            THREE_VECTOR FdTdRk = bioScale(tempK, w);
            THREE_VECTOR FdTdRl = bioScale(tempL, w);


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

            e->virial.xx += virial.xx*ktors;
            e->virial.xy += virial.xy*ktors;
            e->virial.xz += virial.xz*ktors;
            e->virial.yy += virial.yy*ktors;
            e->virial.yz += virial.yz*ktors;
            e->virial.zz += virial.zz*ktors;
        }
    }
    return etorstot;
}

double resImproperSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, ETYPE *e, BIOWEIGHTS* weights)
{
    double eimprtot = 0;

    IMPR_CONN ** imprList = resiConn->imprList;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE imprRange = resiConn->atmRanges[i]->imprRange;
        if (imprRange.start == -1)
        {
            continue;
        }
        for (int imprIndex = imprRange.start; imprIndex < imprRange.end; imprIndex++)
        {
            IMPR_CONN* imprConn = imprList[imprIndex];
            int indexI = resRange->start + imprConn->atmI;
            int indexJ = resRange->start + imprConn->atmJ;
            int indexK = resRange->start + imprConn->atmK;
            int indexL = resRange->start + imprConn->atmL;

            //if not local atom than exit out
            if (gidOrder[indexJ].id >= (int) nlocal)
            {
                continue;
            }
            //if (( state->rx[indexI]==-1 || state->rx[indexJ]==-1 || state->rx[indexK]==-1 || state->rx[indexL]==-1)){continue;}    //FIX ME! 
            //Energy
            //imporper++;

            double impr;
            double sinX;
            //THREE_VECTOR v_ij;
            //THREE_VECTOR v_kj;
            //THREE_VECTOR v_kl;

            THREE_VECTOR dTdrI;
            THREE_VECTOR dTdrJ;
            THREE_VECTOR dTdrK;
            THREE_VECTOR dTdrL;

            THREE_SMATRIX virial;

            if (bioDihedralFast(state, indexI, indexJ, indexK, indexL, &impr, &sinX,
                                &dTdrI, &dTdrJ, &dTdrK, &dTdrL, &virial) < 0)
            {
                continue; // skip the torsion with angle approach to 0.
                // also skip the force calculation because a & b approach to 0;
            }

            IMPROPER_PARMS* imprPtr = imprConn->imprPtr;
            //double imprDelta=impr-imprPtr->psi0/RAD2DEG;
            double imprDelta = impr - imprPtr->psi0;
            double eimpr = imprPtr->kpsi * imprDelta*imprDelta;
            //double w = get4Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ], resiConn->atomList[indexK], resiConn->atomList[indexL]);

            double w = get4Weights(weights, resiConn->atomList[imprConn->atmI], resiConn->atomList[imprConn->atmJ], resiConn->atomList[imprConn->atmK], resiConn->atomList[imprConn->atmL]);
            eimpr *= w;
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
            //kimpr=-kimpr; // Flip the sign for bioVec
            kimpr *= w;
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

            e->virial.xx += virial.xx*kimpr;
            e->virial.xy += virial.xy*kimpr;
            e->virial.xz += virial.xz*kimpr;
            e->virial.yy += virial.yy*kimpr;
            e->virial.yz += virial.yz*kimpr;
            e->virial.zz += virial.zz*kimpr;

        }
    }
    return eimprtot;
}

double resBpairSortedWeighted(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, RESRANGE* resRange, RESI_CONN* resiConn, CHARMMPOT_PARMS *parms, ETYPE *e, double *eblj, double *ebele, BIOWEIGHTS* weights)
{
    *eblj = 0;
    *ebele = 0;
    double ebpairtot = 0;
    double rcut = parms->rmax;
    BPAIR_CONN ** bpairList = resiConn->bpairList;
    for (int i = 0; i < resiConn->atomListSize; i++)
    {
        //for this atom, get the range of adjacent bonds where this atom is listed first in pair data structure
        RANGE bpairRange = resiConn->atmRanges[i]->bpairRange;
        if (bpairRange.start == -1)
        {
            continue;
        } // -1 means no relevant bonds...so go to next pair
        //bondRanges are the intra-residue offsets of each atom in the bond pair...use to get the global indeces 
        for (int bpairIndex = bpairRange.start; bpairIndex < bpairRange.end; bpairIndex++)
        {
            BPAIR_CONN* bpairConn = bpairList[bpairIndex];
            int indexI = resRange->start + (bpairConn->atmI);
            int indexJ = resRange->start + (bpairConn->atmJ);
            if (gidOrder[indexI].id >= (int) nlocal)
            {
                continue;
            }

            // Energy
            // CHARMM LJ: V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
            //double r=bioBond(state, indexI, indexJ);
            THREE_VECTOR rVec;
            bioVec(state, indexI, indexJ, &rVec);
            double r = bioNorm(rVec);

            if (r > rcut)
            {
                continue; // if the pair is larger than the cutoff, don't count it;
            }

            if (r < 0.0001)
            {
                printf("Bond Pair is too close %s-%s, %lu-%lu, %d-%d : %f",
                       state->species[indexI]->name, state->species[indexJ]->name,
                       state->label[indexI], state->label[indexJ],
                       indexI, indexJ, r);
            }
            // L-J term
            double sigma = bpairConn->sigma;
            double eps = bpairConn->eps;
            double ir = 1 / r;
            double sigma_r = sigma * ir;
            double s2 = sigma_r*sigma_r;
            double s4 = s2*s2;
            double s6 = s4*s2;
            double s12 = s6*s6;
            double dvdr = 12.0 * eps * (s6 - s12) * ir;
            double ebpair = eps * (s12 - 2.0 * s6);
            //double w = get2Weights(weights, resiConn->atomList[indexI], resiConn->atomList[indexJ]);

            double w = get2Weights(weights, resiConn->atomList[bpairConn->atmI], resiConn->atomList[bpairConn->atmJ]);
            ebpair *= w;
            ebpairtot = ebpairtot + ebpair;
            *eblj = *eblj + ebpair;

            // Electrostatic term.
            double qI = ((ATOMTYPE_PARMS*) (state->species[indexI])->parm)->charge;
            double qJ = ((ATOMTYPE_PARMS*) (state->species[indexJ])->parm)->charge;
            double ebelec = ke * qI * qJ*ir;
            double dedr = -ebelec*ir; // -ke*qI*qJ/r^2

            ebelec *= w;
            ebpairtot = ebpairtot + ebelec;
            *ebele = *ebele + ebelec;
            /*
                    if(ebpair>1.0){
                        printf("I=%d J=%d SpecieI=%s SpecieJ=%s GidI=%llu, Gid=%llu, BondI=%s BondJ=%s\n", 
                                indexI, indexJ, state->species[indexI]->name, state->species[indexJ]->name,
                                state->label[indexI], state->label[indexJ],
                                bondPtr->atmTypeI, bondPtr->atmTypeJ);         
                        printf("Bond length=%f bondPara=%f bond energy=%f\n",b, bondPtr->b0, eb);
                        printf("x1=%f    y1=%f    z1=%f    x2=%f    y2=%f    z2=%f\n", 
                                state->rx[indexI],state->ry[indexI],state->rz[indexI],
                                state->rx[indexJ],state->ry[indexJ],state->rz[indexJ]);
                    }
             */
            // Force
            THREE_VECTOR unit_ij;
            //bioUnit(state, indexI, indexJ, &unit_ij);
            unit_ij.x = rVec.x*ir;
            unit_ij.y = rVec.y*ir;
            unit_ij.z = rVec.z*ir;

            double kforce = -1 * (dvdr + dedr);
            double fxDelta = w * kforce * unit_ij.x;
            double fyDelta = w * kforce * unit_ij.y;
            double fzDelta = w * kforce * unit_ij.z;
            state->fx[indexI] -= fxDelta; // Fj=-Fi
            state->fy[indexI] -= fyDelta;
            state->fz[indexI] -= fzDelta;
            state->fx[indexJ] += fxDelta; //  F=-dU/dx
            state->fy[indexJ] += fyDelta;
            state->fz[indexJ] += fzDelta;

            // Pressure

            double fpxx = fxDelta * rVec.x;
            double fpxy = fxDelta * rVec.y;
            double fpxz = fxDelta * rVec.z;
            double fpyy = fyDelta * rVec.y;
            double fpyz = fyDelta * rVec.z;
            double fpzz = fzDelta * rVec.z;

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
    }
    //correct for morphing atom charges
    if (resiConn->weightListSize > 0)
    {
        for (int i = 0; i < resiConn->weightListSize; i++)
        {
            int indexI = resRange->start + resiConn->weightList[i]->atmI;
            int indexJ = resRange->start + resiConn->weightList[i]->atmJ;
            THREE_VECTOR rVec;
            bioVec(state, indexI, indexJ, &rVec);
            double r = bioNorm(rVec);
            if (r <= rcut)
            {
                if (r < 0.0001)
                {
                    printf("Non-interacting Pair is too close %s-%s, %lu-%lu, %d-%d : %f",
                           state->species[indexI]->name, state->species[indexJ]->name,
                           state->label[indexI], state->label[indexJ],
                           indexI, indexJ, r);
                }
                double ir = 1 / r;

                //charges should already be weighted 
                // Electrostatic term.
                double qI = ((ATOMTYPE_PARMS*) (state->species[indexI])->parm)->charge;
                double qJ = ((ATOMTYPE_PARMS*) (state->species[indexJ])->parm)->charge;
                double ebelec = ke * qI * qJ*ir;
                double dedr = -ebelec*ir; // -ke*qI*qJ/r^2

                ebpairtot -= ebelec;
                *ebele = *ebele - ebelec;
                // Force
                THREE_VECTOR unit_ij;
                //bioUnit(state, indexI, indexJ, &unit_ij);
                unit_ij.x = rVec.x*ir;
                unit_ij.y = rVec.y*ir;
                unit_ij.z = rVec.z*ir;

                double kforce = -1 * dedr;
                double fxDelta = kforce * unit_ij.x;
                double fyDelta = kforce * unit_ij.y;
                double fzDelta = kforce * unit_ij.z;
                state->fx[indexI] += fxDelta; // Fj=-Fi
                state->fy[indexI] += fyDelta;
                state->fz[indexI] += fzDelta;
                state->fx[indexJ] -= fxDelta; //  F=-dU/dx
                state->fy[indexJ] -= fyDelta;
                state->fz[indexJ] -= fzDelta;

                // Pressure

                double fpxx = fxDelta * rVec.x;
                double fpxy = fxDelta * rVec.y;
                double fpxz = fxDelta * rVec.z;
                double fpyy = fyDelta * rVec.y;
                double fpyz = fyDelta * rVec.z;
                double fpzz = fzDelta * rVec.z;

                e->virial.xx -= fpxx; // appears to be a sign error. 
                e->virial.xy -= fpxy;
                e->virial.xz -= fpxz;
                e->virial.yy -= fpyy;
                e->virial.yz -= fpyz;
                e->virial.zz -= fpzz;
            }
        }
    }

    return ebpairtot;
}

