#include "nglfrattle.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "box.h"
#include "species.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "expandbuffer.h"
#include "auxNeighbor.h"
#include "preduce.h"
#include "group.h"
#include "state.h"

#include "bioCharmm.h"
#include "bioCharmmParms.h"
//#include "setList.h"

void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
void nglfrattle_collision(SYSTEM *sys, double dt);
double pair_function(SYSTEM *sys, int i, int j, double r, double *dr);


static AUXNEIGHBOR *auxNeighbor;

NGLFRATTLE_PARMS *nglfrattle_parms(INTEGRATOR*integrator)
{
    NGLFRATTLE_PARMS *parms;
    parms = ddcMalloc(sizeof (NGLFRATTLE_PARMS));
    parms->acc = NULL;
    double rmax = 3.0;
    auxNeighbor = auxNeighbor_request(rmax);
    // Get charmmpotparms initialized
    SIMULATE *simulate = (SIMULATE*) integrator->parent;
    SYSTEM *system = simulate->system;
    int foundCHARMM = 0;
    for (int i = 0; i < system->npotential; ++i)
    {
        //if (strcmp(system->potential[i]->type, "CHARMM") == 0){
        if (system->potential[i]->itype == CHARMM)
        {
            parms->charmmpot_parms = (CHARMMPOT_PARMS *) (system->potential[i]->parms);
            foundCHARMM = 1;
            break;
        }
    }
    if (foundCHARMM < 0)
    {
        printf("nglfrattle_parms: Cannot find CHARMM parameter for nglfrattle\n");
    }
    return parms;
}

void scalePositionsByBoxChange_sk(BOX_STRUCT *box, double time, double *rx, double *ry, double *rz, unsigned nlocal)
{
    THREE_MATRIX hfac;
    THREE_VECTOR rold, r;
    box_put(box, BOX_TIME, (void *) &time);
    box_get(box, HFAC, (void *) &hfac);
    if (!matrix_equal(hfac, I_3x3))
    {
        for (unsigned kk = 0; kk < nlocal; kk++)
        {
            rold.x = rx[kk];
            rold.y = ry[kk];
            rold.z = rz[kk];
            r = matrix_vector(hfac, rold);
            rx[kk] = r.x;
            ry[kk] = r.y;
            rz[kk] = r.z;

        }
    }
}

void updateStateAliases_sk(SYSTEM *sys, unsigned *nlocal, unsigned *nion, double **rx, double **ry, double **rz, double **vx, double **vy, double **vz, double **fx, double **fy, double **fz, SPECIES ***species, gid_type **label)
{
    STATE *state = sys->collection->state;
    *nlocal = sys->nlocal;
    *nion = sys->nion;
    *rx = state->rx;
    *ry = state->ry;
    *rz = state->rz; // The SYSTEM and STATE might change during the call to ddcenergy
    *vx = state->vx;
    *vy = state->vy;
    *vz = state->vz; // (i.e., we might reassign particles to new tasks) so we need to
    *fx = state->fx;
    *fy = state->fy;
    *fz = state->fz; // update all of the aliases we use.
    *label = state->label;
    *species = state->species;
}

int compareSoloGid(const void *v1, const void *v2)
{ //Name of compareGID is taken.
    const gid_type *u1 = (gid_type *) v1;
    const gid_type *u2 = (gid_type *) v2;
    if (*u1<*u2)
    {
        return -1;
    }
    else if (*u1>*u2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void resMoveA(SIMULATE*simulate, STATE* state, GID_ORDER* gidOrder, int nTotal,
              char * name, RESRANGE resRange, RESI_CONN* resiConn)
{


    const int debug = 1;

    double dt = simulate->dt;
    //double dt2=dt/2.0;
    //double dtsq2=dt*dt2;
    //double time = simulate->time;

    int nlocal = simulate->system->nlocal;

    const double RPTOL = 1.0e-6;
    const double tol = 1.0e-5; // unit less
    const double tol2 = 2 * tol;
    const int maxit = 500;

    //RESI_CONN* resiConn=findResiConnBySpecieName(charmmParms, name);

    REF *v0 = vaf_v0();

    for (int i = 0; i < resiConn->heavyAtmsSize; i++)
    {
        HEAVYATM* heavyAtm = resiConn->heavyAtms[i];
        int indexHeavyAtm = resRange.start + heavyAtm->atomID;
        if (indexHeavyAtm > nTotal)
        {
            continue;
        } //skip if atom doesn't exist
        if (state->rx[indexHeavyAtm] == -1)
        {
            continue;
        } //FIX ME!
        if (gidOrder[indexHeavyAtm].id >= (int) nlocal)
        {
            continue;
        } //Not in local

        if (debug)
        {
            //double check if the bonded hydrogen atoms is in the local+remote
            for (int j = 0; j < heavyAtm->numHydrogen; j++)
            {
                int indexH = resRange.start + heavyAtm->hydrogen[j];
                if (indexH > nTotal)
                {
                    printf("resMove: hydrogen doesn't exist.\n");
                }
                if (state->rx[indexH] == -1)
                {
                    printf("resMove: hydrogen .\n"); // Ask Shiv
                }
            }
        }

        if (heavyAtm->numHydrogen == 0)
        {
            // Velocity Verlet without constraints

            THREE_VECTOR delta;
            state->rx[indexHeavyAtm] += delta.x = dt * state->vx[indexHeavyAtm];
            state->ry[indexHeavyAtm] += delta.y = dt * state->vy[indexHeavyAtm];
            state->rz[indexHeavyAtm] += delta.z = dt * state->vz[indexHeavyAtm];
            if (v0 != NULL) VOP1(v0[indexHeavyAtm].r, +=, delta);
        }
        else
        {
            // Velocity Verlet with constraints
            int numAtm = heavyAtm->numHydrogen + 1; // total number of atom 

            double *rx = ddcMalloc(numAtm * sizeof (double));
            double *ry = ddcMalloc(numAtm * sizeof (double));
            double *rz = ddcMalloc(numAtm * sizeof (double));

            double *vx = ddcMalloc(numAtm * sizeof (double));
            double *vy = ddcMalloc(numAtm * sizeof (double));
            double *vz = ddcMalloc(numAtm * sizeof (double));

            double *px = ddcMalloc(numAtm * sizeof (double));
            double *py = ddcMalloc(numAtm * sizeof (double));
            double *pz = ddcMalloc(numAtm * sizeof (double));

            double *mass = ddcMalloc(numAtm * sizeof (double));
            double *bond2 = ddcMalloc(numAtm * sizeof (double));

            for (int j = 0; j < numAtm; ++j)
            {
                int index = 0;
                if (j == 0)
                {
                    // First elements of arrays are for heavy atom;
                    index = indexHeavyAtm;
                    bond2[j] = 0; //leave the first element empty
                }
                else
                {
                    // The rest elements are hydrogen
                    index = resRange.start + heavyAtm->hydrogen[j - 1];
                    double bond = heavyAtm->hBndConn[j - 1]->bondPtr->b0;
                    bond2[j] = bond*bond;
                }
                mass[j] = ((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
                //double ax=state->fx[index]/mass[j];
                //double ay=state->fy[index]/mass[j];
                //double az=state->fz[index]/mass[j];
                rx[j] = state->rx[index];
                ry[j] = state->ry[index];
                rz[j] = state->rz[index];
                vx[j] = state->vx[index];
                vy[j] = state->vy[index];
                vz[j] = state->vz[index];
                px[j] = rx[j] + vx[j] * dt;
                py[j] = ry[j] + vy[j] * dt;
                pz[j] = rz[j] + vz[j] * dt;

            }

            int *moving = ddcMalloc(numAtm * sizeof (int));
            int *moved = ddcMalloc(numAtm * sizeof (int));
            for (int j = 0; j < numAtm; ++j)
            {
                moving[j] = 0;
                moved[j] = 1;
            }

            int it = 0;
            int DONE = 0;
            //start the iterative loop
            while (DONE == 0 && it < maxit)
            {
                DONE = 1;

                for (int j = 1; j < numAtm; ++j)
                {
                    if (moved[0] == 1 || moved[j] == 1)
                    {
                        double pxab = px[0] - px[j];
                        double pyab = py[0] - py[j];
                        double pzab = pz[0] - pz[j];
                        nearestImage(&pxab, &pyab, &pzab);
                        double pab2 = pxab * pxab + pyab * pyab + pzab*pzab;
                        double diffab2 = bond2[j] - pab2;
                        if (fabs(diffab2) > bond2[j] * tol2)
                        {
                            double rxab = rx[0] - rx[j];
                            double ryab = ry[0] - ry[j];
                            double rzab = rz[0] - rz[j];
                            nearestImage(&rxab, &ryab, &rzab);
                            double rpab = rxab * pxab + ryab * pyab + rzab*pzab;
                            if (rpab < bond2[j] * RPTOL)
                            {
                                printf("resMoveA: Constraint failure\n");
                                continue;
                            }

                            double rma = 1.0 / mass[0];
                            double rmb = 1.0 / mass[j];
                            double gab = diffab2 / (2.0 * (rma + rmb) * rpab);
                            double dx = rxab*gab;
                            double dy = ryab*gab;
                            double dz = rzab*gab;

                            //printf("J=%d\n", j);
                            //printf("Before: px[0]=%f py[0]=%f pz[0]=%f\n", px[0], py[0], pz[0]);
                            //printf("Before: px[j]=%f py[j]=%f pz[j]=%f\n", px[j], py[j], pz[j]);
                            px[0] = px[0] + rma*dx;
                            py[0] = py[0] + rma*dy;
                            pz[0] = pz[0] + rma*dz;
                            px[j] = px[j] - rmb*dx;
                            py[j] = py[j] - rmb*dy;
                            pz[j] = pz[j] - rmb*dz;

                            //printf("After: px[0]=%f py[0]=%f pz[0]=%f\n", px[0], py[0], pz[0]);
                            //printf("After: px[j]=%f py[j]=%f pz[j]=%f\n", px[j], py[j], pz[j]);                            
                            //printf("During dx=%f dy=%f dz=%f\n", dx, dy, dz);

                            double vdx = dx / dt;
                            double vdy = dy / dt;
                            double vdz = dz / dt;

                            vx[0] = vx[0] + rma*vdx;
                            vy[0] = vy[0] + rma*vdy;
                            vz[0] = vz[0] + rma*vdz;
                            vx[j] = vx[j] - rmb*vdx;
                            vy[j] = vy[j] - rmb*vdy;
                            vz[j] = vz[j] - rmb*vdz;

                            moving[0] = 1;
                            moving[j] = 1;
                            DONE = 0;
                        }
                    }
                }

                for (int j = 0; j < numAtm; ++j)
                {
                    moved[j] = moving[j];
                    moving[j] = 0;
                }

                it = it + 1;
            }
            //end of iterative loop
            //printf("resname=%s heavy atom=%d iteration=%d\n", name, i, it);

            if (DONE == 0)
            {
                printf("resMoveA: too many constraint iterations.\n");
            }

            // Store the new values
            for (int j = 0; j < numAtm; ++j)
            {
                int index = 0;
                if (j == 0)
                {
                    // First elements of arrays are for heavy atom;
                    index = indexHeavyAtm;
                }
                else
                {
                    // The rest elements are hydrogen
                    index = resRange.start + heavyAtm->hydrogen[j - 1];
                }

                if (v0 != NULL)
                {
                    THREE_VECTOR delta;
                    delta.x = px[j] - state->rx[index];
                    delta.y = py[j] - state->ry[index];
                    delta.z = pz[j] - state->rz[index];
                    VOP1(v0[index].r, +=, delta);
                }

                state->rx[index] = px[j];
                state->ry[index] = py[j];
                state->rz[index] = pz[j];

                state->vx[index] = vx[j];
                state->vy[index] = vy[j];
                state->vz[index] = vz[j];
            }

            ddcFree(rx);
            ddcFree(ry);
            ddcFree(rz);
            ddcFree(vx);
            ddcFree(vy);
            ddcFree(vz);
            ddcFree(px);
            ddcFree(py);
            ddcFree(pz);
            ddcFree(mass);
            ddcFree(bond2);
            ddcFree(moving);
            ddcFree(moved);
        }
    }

}

void resMoveB(SIMULATE*simulate, STATE* state, GID_ORDER* gidOrder, int nTotal,
              char * name, RESRANGE resRange, RESI_CONN* resiConn)
{


    const int debug = 1;

    double dt = simulate->dt;
    double dt2 = dt / 2.0;

    double k = 0.0;
    double wc = 0.0;

    //double time = simulate->time;
    int nlocal = simulate->system->nlocal;

    const double tol = 1.0e-5; //unit less
    const int maxit = 500;

    //RESI_CONN* resiConn=findResiConnBySpecieName(charmmParms, name);

    for (int i = 0; i < resiConn->heavyAtmsSize; i++)
    {
        HEAVYATM* heavyAtm = resiConn->heavyAtms[i];
        int indexHeavyAtm = resRange.start + heavyAtm->atomID;
        if (indexHeavyAtm > nTotal)
        {
            continue;
        } //skip if atom doesn't exist
        if (state->rx[indexHeavyAtm] == -1)
        {
            continue;
        } //FIX ME!
        if (gidOrder[indexHeavyAtm].id >= (int) nlocal)
        {
            continue;
        } //Not in local

        if (debug)
        {
            //double check if the bonded hydrogen atoms is in the local+remote
            for (int j = 0; j < heavyAtm->numHydrogen; j++)
            {
                int indexH = resRange.start + heavyAtm->hydrogen[j];
                if (indexH > nTotal)
                {
                    printf("resMove: hydrogen doesn't exist.\n");
                }
                if (state->rx[indexH] == -1)
                {
                    printf("resMove: hydrogen .\n"); // Ask Shiv
                }
            }
        }

        if (heavyAtm->numHydrogen != 0)
        {
            // Velocity Verlet with constraints
            int numAtm = heavyAtm->numHydrogen + 1; // total number of atom 

            double *rx = ddcMalloc(numAtm * sizeof (double));
            double *ry = ddcMalloc(numAtm * sizeof (double));
            double *rz = ddcMalloc(numAtm * sizeof (double));

            double *vx = ddcMalloc(numAtm * sizeof (double));
            double *vy = ddcMalloc(numAtm * sizeof (double));
            double *vz = ddcMalloc(numAtm * sizeof (double));

            double *mass = ddcMalloc(numAtm * sizeof (double));

            double * bond2 = ddcMalloc(numAtm * sizeof (double));


            for (int j = 0; j < numAtm; ++j)
            {
                int index = 0;
                if (j == 0)
                {
                    // First elements of arrays are for heavy atom;
                    index = indexHeavyAtm;
                    bond2[j] = 0; //leave the first element empty
                }
                else
                {
                    // The rest elements are hydrogen
                    index = resRange.start + heavyAtm->hydrogen[j - 1];
                    double bond = heavyAtm->hBndConn[j - 1]->bondPtr->b0;
                    bond2[j] = bond*bond;
                }
                mass[j] = ((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
                rx[j] = state->rx[index];
                ry[j] = state->ry[index];
                rz[j] = state->rz[index];
                vx[j] = state->vx[index];
                vy[j] = state->vy[index];
                vz[j] = state->vz[index];
            }

            int *moving = ddcMalloc(numAtm * sizeof (int));
            int *moved = ddcMalloc(numAtm * sizeof (int));
            for (int j = 0; j < numAtm; ++j)
            {
                moving[j] = 0;
                moved[j] = 1;
            }

            int it = 0;
            int DONE = 0;
            //start the iterative loop
            while (DONE == 0 && it < maxit)
            {
                DONE = 1;

                for (int j = 1; j < numAtm; ++j)
                {
                    if (moved[0] == 1 || moved[j] == 1)
                    {
                        double vxab = vx[0] - vx[j];
                        double vyab = vy[0] - vy[j];
                        double vzab = vz[0] - vz[j];

                        double rxab = rx[0] - rx[j];
                        double ryab = ry[0] - ry[j];
                        double rzab = rz[0] - rz[j];

                        nearestImage(&rxab, &ryab, &rzab);

                        double rvab = rxab * vxab + ryab * vyab + rzab*vzab;
                        double rma = 1.0 / mass[0];
                        double rmb = 1.0 / mass[j];
                        double gab = -rvab / ((rma + rmb) * bond2[j]);

                        if (fabs(gab) > tol)
                        {
                            wc = wc + gab * bond2[j];
                            double dx = rxab*gab;
                            double dy = ryab*gab;
                            double dz = rzab*gab;

                            vx[0] = vx[0] + rma*dx;
                            vy[0] = vy[0] + rma*dy;
                            vz[0] = vz[0] + rma*dz;
                            vx[j] = vx[j] - rmb*dx;
                            vy[j] = vy[j] - rmb*dy;
                            vz[j] = vz[j] - rmb*dz;

                            moving[0] = 1;
                            moving[j] = 1;
                            DONE = 0;
                        }
                    }
                }

                for (int j = 0; j < numAtm; ++j)
                {
                    moved[j] = moving[j];
                    moving[j] = 0;
                }

                it = it + 1;
            }
            //end of iterative loop

            if (DONE == 0)
            {
                printf("resMoveB: too many contraint iterations.\n");
            }

            // Store the new values
            for (int j = 0; j < numAtm; ++j)
            {
                int index = 0;
                if (j == 0)
                {
                    // First elements of arrays are for heavy atom;
                    index = indexHeavyAtm;
                }
                else
                {
                    // The rest elements are hydrogen
                    index = resRange.start + heavyAtm->hydrogen[j - 1];
                }
                state->vx[index] = vx[j];
                state->vy[index] = vy[j];
                state->vz[index] = vz[j];

                k = k + mass[j]*(vx[j] * vx[j] + vy[j] * vy[j] + vz[j] * vz[j]);
            }

            ddcFree(rx);
            ddcFree(ry);
            ddcFree(rz);
            ddcFree(vx);
            ddcFree(vy);
            ddcFree(vz);
            ddcFree(mass);
            ddcFree(bond2);
            ddcFree(moving);
            ddcFree(moved);
        }

    }

    k = k * 0.5; // TODO: do we need k here?
    wc = wc / dt2 / 3.0; // TODO: do we need wc here?

}

void padding(DDC*ddc, SIMULATE*simulate, NGLFRATTLE_PARMS*p)
{

    SYSTEM* sys = simulate->system;
    STATE* state = sys->collection->state;
    CHARMMPOT_PARMS *charmmpot_parms = p->charmmpot_parms;
    //CHARMM_PARMS *charmmParms=charmmpot_parms->charmmParms;

    SETLIST *residueSet = &charmmpot_parms->residueSet;
    LISTNODE* residueList = residueSet->list;

    GID_ORDER gidOrder2[residueSet->molSize];
    GID_ORDER *gidOrder = charmmpot_parms->gidOrder;
    unsigned partialNion = residueSet->molSize;
    int exists[residueSet->molSize];

    STATE *statechpad = &(charmmpot_parms->statechpad);
    statechpad->nlocal = state->nlocal;
    statechpad->nion = partialNion;

    unsigned nion = sys->nion;

    unsigned s_id = 0;
    int sp_id = 0;
    for (int i = nion; i < residueSet->molSize; i++)
    {
        gidOrder2[i].id = INT_MAX; //this will be removed in next commit. see note below on why we set gid to this high number. 
        statechpad->label[i] = INT_MAX;
        ;
    }
    unsigned index = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN * resiConn = residueList[i].resiConn;
        char * name = resiConn->resName;
        char *terminus = "x";
        if (resiConn->nTer)
        {
            terminus = "n";
        }
        if (resiConn->cTer)
        {
            terminus = "c";
        }
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            //char atm[11];
            //char *atmName = atm;
            char atmName[11];
            strcpy(atmName, name);
            strcat(atmName, terminus);
            strcat(atmName, resiConn->atomList[j]->atmName);
            statechpad->rx[sp_id] = -1;
            s_id = gidOrder[index].id;
            if (index < sys->nion && state->species[s_id] && (state->species[s_id]->name) &&(strcmp(atmName, state->species[s_id]->name) == 0))
            {
                if (statechpad->label[sp_id] != state->label[s_id])
                {
                    printf("pad gid=%lu  reg gid=%lu\n", statechpad->label[sp_id], state->label[s_id]);
                }
                statechpad->rx[sp_id] = state->rx[s_id];
                statechpad->ry[sp_id] = state->ry[s_id];
                statechpad->rz[sp_id] = state->rz[s_id];
                statechpad->vx[sp_id] = state->vx[s_id];
                statechpad->vy[sp_id] = state->vy[s_id];
                statechpad->vz[sp_id] = state->vz[s_id];
                statechpad->fx[sp_id] = state->fx[s_id];
                statechpad->fy[sp_id] = state->fy[s_id];
                statechpad->fz[sp_id] = state->fz[s_id];
                statechpad->species[sp_id] = state->species[s_id];
                gidOrder2[sp_id].id = gidOrder[index].id;
                exists[sp_id] = 1;
                index++;
            }
            else
            { //we have a hole in the padded array, some residue atoms are missing
                //set gidOrder2 to
                //unrealistically high value: we use this info in energy functions so we know when to ignore terms
                //but we can just use the exists array mentioned below
                //this will be removed in next commit
                statechpad->label[sp_id] = INT_MAX;
                gidOrder2[sp_id].id = INT_MAX;
                exists[sp_id] = 0;
            }
            sp_id++;
        }
    }

    //RESRANGE resRangeTmp[count];
    RESRANGE resRange2[residueSet->listSize];
    unsigned start = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        resRange2[i].start = start;
        start += residueList[i].size;
        resRange2[i].end = start;
    }

    //CHARMM_PARMS* charmmParms = p->charmmpot_parms->charmmParms;

    for (int i = 0; i < residueSet->listSize; i++)
    {
        p->moveFunc(simulate, statechpad, gidOrder2, residueSet->molSize,
                    residueList[i].name, resRange2[i], residueList[i].resiConn);

    }

    // Update the force to ddc items 
    // May not need to update coordinates.
    // transfer state from padded sorted state array to system state array
    for (int i = 0; i < residueSet->molSize; i++)
    {
        if (!exists[i])
        {
            continue;
        };
        unsigned index = gidOrder2[i].id;
        if (index > sys->nion)
        {
            continue;
        }
        state->rx[index] = statechpad->rx[i];
        state->ry[index] = statechpad->ry[i];
        state->rz[index] = statechpad->rz[i];
        state->vx[index] = statechpad->vx[i];
        state->vy[index] = statechpad->vy[i];
        state->vz[index] = statechpad->vz[i];
        state->fx[index] = statechpad->fx[i];
        state->fy[index] = statechpad->fy[i];
        state->fz[index] = statechpad->fz[i];
    }

    profile(CHARMM_COVALENT, END);

}

void nglfrattle(DDC*ddc, SIMULATE*simulate, NGLFRATTLE_PARMS*p)
{
    double dt = simulate->dt;

    double time = simulate->time;
    SYSTEM* sys = simulate->system;
    STATE* state = sys->collection->state;
    //CHARMMPOT_PARMS *charmmpot_parms = p->charmmpot_parms;
    unsigned nlocal;
    unsigned nion;
    double *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz;
    SPECIES** species;
    LONG64* label;

    updateStateAliases_sk(sys, &nlocal, &nion, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    scalePositionsByBoxChange_sk(sys->box, time, rx, ry, rz, nlocal);

    for (unsigned kk = 0; kk < nlocal; kk++)
    {
        GROUP* group = state->group[kk];
        group->velocityUpdate(FRONT_TIMESTEP, kk, group, state, time, 0.5 * dt);
    }

    p->moveFunc = (void(*)(SIMULATE*, STATE*, GID_ORDER*, int, char *, RESRANGE, CHARMM_PARMS*))resMoveA;
    padding(ddc, simulate, p);

    double time_plus_dt = time + dt;
    scalePositionsByBoxChange_sk(sys->box, time_plus_dt, rx, ry, rz, nlocal);
    for (unsigned kk = 0; kk < nlocal; kk++)
    {
        backInBox_fast(rx + kk, ry + kk, rz + kk);
    }

    ddc->update = 0;
    time += dt; // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
    simulate->time = sys->time = time;
    simulate->loop++;
    sys->loop = simulate->loop;

    //group_Update function is an empty function for free .

    for (int kk = 0; kk < sys->ngroup; kk++)
    {
        sys->group[kk]->Update1(sys->group[kk], -1, state, time, 0.5 * dt);
    }

    if (ddcenergy(ddc, sys, 0) != 0) return;
    updateStateAliases_sk(sys, &nlocal, &nion, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    for (int kk = 0; kk < sys->ngroup; kk++)
    {
        sys->group[kk]->Update(sys->group[kk], BACK_TIMESTEP, state, time, 0.5 * dt);
    }

    for (unsigned kk = 0; kk < nlocal; kk++)
    {
        GROUP* group = state->group[kk];
        group->velocityUpdate(BACK_TIMESTEP, kk, group, state, time, 0.50 * dt);
    }

    p->moveFunc = (void(*)(SIMULATE*, STATE*, GID_ORDER*, int, char *, RESRANGE, CHARMM_PARMS*))resMoveB;
    padding(ddc, simulate, p);

    kinetic_terms(sys, 1);

    for (int kk = 0; kk < sys->ngroup; kk++)
    {
        sys->group[kk]->Update2(sys->group[kk], -1, state, time, 0.5 * dt);
    }
    //eval_energyInfo(sys);
    for (int kk = 0; kk < sys->ngroup; kk++)
    {
        sys->group[kk]->Update(sys->group[kk], FRONT_TIMESTEP, state, time, 0.5 * dt);
    }

    /*errorCheck(ddc->domain_id, simulate->loop, state, sys->energyInfo, p, datafile); */
    simulate->time = sys->time = time;

}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
