#include "nglfconstraint.h"
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
#include "mpiUtils.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "codata.h"
#include "units.h"

void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
void nglfconstraint_collision(SYSTEM *sys, double dt);
double pair_function(SYSTEM *sys, int i, int j, double r, double *dr);
char *strcat_(char *a, char *b);
void velocityConstraint(double dt, STATE *state, NGLFCONSTRAINT_PARMS *p, GID_ORDER *gidOrder2, RESRANGE *resRange2, int);
double frontFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab);
double backFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab);
void solve(int n, double a[n][n], double *b, double *x);
void solveTest(int n, double a[n][n], double *b, double *x);
void copyState2Statechpad(STATE *state, STATE *statechpad, GID_ORDER *gidOrder2, SETLIST *residueSet);
//THREE_SMATRIX molecularVirial(SYSTEM *sys,THREE_SMATRIX virial);
THREE_SMATRIX molecularPressure(SYSTEM *sys, THREE_SMATRIX virial, double T);

static AUXNEIGHBOR *auxNeighbor;

void adjustPosn(STATE *state, BOX_STRUCT* box)
{
    THREE_MATRIX hfac;
    box_get(box, HFAC, (void *) &hfac);
    for (int kk = 0; kk < state->nlocal; kk++)
    {
        THREE_VECTOR old = {state->rx[kk], state->ry[kk], state->rz[kk]};
        THREE_VECTOR new = matrix_vector(hfac, old);
        state->rx[kk] = new.x;
        state->ry[kk] = new.y;
        state->rz[kk] = new.z;
    }
}


/* Local Variables: */
/* tab-width: 3 */

/* End: */
void changeVolume(STATE *state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt)
{
    double btt = beta * dt / tau;
    double Pxx, Pyy, Pzz;
    Pxx = Pyy = 0.5 * (pTensor->xx + pTensor->yy);
    Pzz = pTensor->zz;
    //Pxx = Pyy = Pzz=(1.0/3.0)*(pressureTensor->xx+pressureTensor->yy+pressureTensor->zz);   // Isotropic; 

    THREE_MATRIX lambda = mzero;
    lambda.xx = cbrt(1.0 + Pxx * btt);
    lambda.yy = cbrt(1.0 + Pyy * btt);
    //lambda.xx = 1.0;
    //lambda.yy = 1.0;
    lambda.zz = cbrt(1.0 + Pzz * btt);
    THREE_MATRIX h0 = box_get_h(box);
    //double time =system_getTime(NULL);
    THREE_MATRIX h = matrix_matrix(lambda, h0);
    box_put(box, HO, &h);
    adjustPosn(state, box);

}

NGLFCONSTRAINT_PARMS *nglfconstraint_parms(INTEGRATOR*integrator)
{
    NGLFCONSTRAINT_PARMS *parms;
    parms = ddcMalloc(sizeof (NGLFCONSTRAINT_PARMS));
    object_get((OBJECT*) integrator, "T", &parms->T, WITH_UNITS, 1, "310", "T", NULL);
    object_get((OBJECT*) integrator, "P0", &parms->P0, WITH_UNITS, 1, "0.0", "pressure", NULL);
    object_get((OBJECT*) integrator, "beta", &parms->beta, WITH_UNITS, 1, "0.0", "1/pressure", NULL);
    object_get((OBJECT*) integrator, "tauBarostat", &parms->tauBarostat, WITH_UNITS, 1, "0.0", "t", NULL);
    object_get((OBJECT*) integrator, "isotropic", &parms->isotropic, INT, 1, "0");
    parms->acc = NULL;
    double rmax = 3.0;
    auxNeighbor = auxNeighbor_request(rmax);
    // Get charmmpotparms initialized
    SIMULATE *simulate = (SIMULATE*) integrator->parent;
    SYSTEM *system = simulate->system;
    int foundMMFF = 0;
    for (int i = 0; i < system->npotential; ++i)
    {
        if (system->potential[i]->itype == CHARMM || system->potential[i]->itype == MARTINI)
        {
            parms->charmmpot_parms = (CHARMMPOT_PARMS *) (system->potential[i]->parms);
            foundMMFF = 1;
            break;
        }
    }
    if (foundMMFF == 0)
        printf("nglfconstraint_parms: Cannot find CHARMM or Martini parameter for nglfconstraint\n");

    if (integrator->uses_gpu == 1)
    {
        nglfconstraintGPU_parms(system, parms);
    }
    return parms;
}
//p2-d2 = (r + v*dt)^2 - d2  = 0.5*(r2-d2)/dt + r.v +0.5*dt*v2 

double frontFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab)
{
    THREE_VECTOR pab = rab;
    VSVOP(pab, +=, dt, *, vab);
    double pab2 = VSQ(pab);
    double rvab = (pab2 - dist2) / (2 * dt);
    //             p2-d2 = (r + v*dt)^2 - d2  = 0.5*(r2-d2)/dt + r.v +0.5*dt*v2 

    return rvab;
}

double backFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab)
{
    double rvab = DOT(rab, vab);
    return rvab;
}

double solveConstraintMatrix(int location, int n, int m, THREE_VECTOR *r, CONS_PAIR **consPair, double *rMass, double dt, double *lambda, THREE_VECTOR *V)
{
    double ( *func) (double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab);
    if (location == FRONT_TIMESTEP) func = frontFunc;
    if (location == BACK_TIMESTEP) func = backFunc;
    double rhs[n], M[n][n];
    for (int ab = 0; ab < n; ab++)
    {
        int a = consPair[ab]->atomIindex;
        int b = consPair[ab]->atomJindex;
        for (int uv = 0; uv < n; uv++)
        {
            int u = consPair[uv]->atomIindex;
            int v = consPair[uv]->atomJindex;
            M[ab][uv] = DOT(r[ab], r[uv])*(((u == a) - (v == a)) * rMass[a] - ((u == b) - (v == b)) * rMass[b]);
        }
        double dist2 = SQ(consPair[ab]->distance);
        THREE_VECTOR Vab;
        VOP2(Vab, =, V[a], -, V[b]);
        double rvab = -func(dt, dist2, r[ab], Vab);
        rhs[ab] = rvab;
    }
    solve(n, M, rhs, lambda);
    for (int a = 0; a < m; a++)
    {
        for (int uv = 0; uv < n; uv++)
        {
            int u = consPair[uv]->atomIindex;
            int v = consPair[uv]->atomJindex;
            double Iabuv = ((u == a) - (v == a)) * rMass[a];
            VSVOP(V[a], +=, Iabuv * lambda[uv], *, r[uv]);
        }
    }
    double err = 0;
    if (location == BACK_TIMESTEP) return err;
    for (int ab = 0; ab < n; ab++) err += (lambda[ab] * lambda[ab]);
    err = sqrt(err);
    return err;
}
//void resMoveCons(int  iResi, int consListSize, CONSTRAINT **consList, double dt, STATE* state, int resRangeStart, int location)

void resMoveConsOld(int consListSize, CONSTRAINT **consList, double dt, STATE* state, THREE_SMATRIX *sion, int resRangeStart, int location)
{

    double ( *func) (double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab);
    if (location == FRONT_TIMESTEP) func = frontFunc;
    if (location == BACK_TIMESTEP) func = backFunc;
    const double tol = 1.0e-12; //unit less
    const int maxit = 500;

    for (int i = 0; i < consListSize; i++)
    {
        CONSTRAINT* constraint = consList[i];
        int numAtm = constraint->numAtom; // total number of atom
        double rMass[numAtm];
        THREE_VECTOR r[numAtm], v[numAtm];
        for (int j = 0; j < numAtm; ++j)
        {
            int index = resRangeStart + constraint->atomIDList[j];
            rMass[j] = 1.0 / ((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
            VSET(r[j], state->rx[index], state->ry[index], state->rz[index]);
            VSET(v[j], state->vx[index], state->vy[index], state->vz[index]);
        }
        int numPair = constraint->numPair;
        THREE_VECTOR rab[numPair];
        double gamma[numPair];
        for (int ab = 0; ab < numPair; ++ab)
        {
            CONS_PAIR* consPair = constraint->conspairList[ab];
            int a = consPair->atomIindex;
            int b = consPair->atomJindex;
            VOP2(rab[ab], =, r[a], -, r[b]);
            nearestImage(&rab[ab].x, &rab[ab].y, &rab[ab].z);
        }
        double errMax = 0.0;
        int it = 0;
        for (int ab = 0; ab < numPair; ++ab) gamma[ab] = 0.0;
        for (; it < maxit; it++) //start the iterative loop
        {
            errMax = 0.0;
            for (int ab = 0; ab < numPair; ++ab)
            {
                CONS_PAIR* consPair = constraint->conspairList[ab];
                int a = consPair->atomIindex;
                int b = consPair->atomJindex;
                double dist2 = SQ(consPair->distance);
                THREE_VECTOR vab;
                VOP2(vab, =, v[a], -, v[b]);
                double rma = rMass[a];
                double rmb = rMass[b];

                double rvab = func(dt, dist2, rab[ab], vab) / dist2;

                double gab = -rvab / (rma + rmb); //units: mass/time)
                double err = fabs(rvab * dt);
                if (err > errMax) errMax = err;
                VSVOP(v[a], +=, (rma * gab), *, rab[ab]);
                VSVOP(v[b], -=, (rmb * gab), *, rab[ab]);
                gamma[ab] += gab;
            }
            if (errMax < tol)
            {
                break;
            }
        } //end of iterative loop
        if (it == maxit) printf("%d resMove: too many contraint iterations.\n", getRank(0));
        *sion = szero;
        for (int ab = 0; ab < numPair; ++ab)
        {
            THREE_VECTOR fab;
            VSVOP(fab, =, (-2.0 / dt * gamma[ab]), *, rab[ab]); //fab = gamma/(dt/2)*rab. delta moment v*m = f *t ;  In this case t = dt/2
            sion->xx += fab.x * rab[ab].x;
            sion->yy += fab.y * rab[ab].y;
            sion->zz += fab.z * rab[ab].z;
            sion->xy += fab.x * rab[ab].y;
            sion->xz += fab.x * rab[ab].z;
            sion->yz += fab.y * rab[ab].z;
        }
        for (int j = 0; j < numAtm; ++j) // Store the new values
        {
            int index = resRangeStart + constraint->atomIDList[j];
            state->vx[index] = v[j].x;
            state->vy[index] = v[j].y;
            state->vz[index] = v[j].z;
        }
    }
}

void resMoveCons(CONSTRAINT *constraint, double dt, STATE* state, int resRangeStart, int location)
{

    const double tol = 1.0e-12; //unit less
    //const int maxit=500;

    //for(int i=0; i<consListSize; i++)
    //{
    // CONSTRAINT* constraint = consList[i]; 
    int numAtm = constraint->numAtom; // total number of atom
    double rMass[numAtm];
    THREE_VECTOR r[numAtm], v[numAtm];
    for (int j = 0; j < numAtm; ++j)
    {
        int index = resRangeStart + constraint->atomIDList[j];
        rMass[j] = 1.0 / ((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
        VSET(r[j], state->rx[index], state->ry[index], state->rz[index]);
        VSET(v[j], state->vx[index], state->vy[index], state->vz[index]);
    }
    int numPair = constraint->numPair;
    //THREE_VECTOR rab[numPair],v0[numPair];
    THREE_VECTOR rab[numPair];
    //double gamma[numPair],lambda[numPair]; 
    double lambda[numPair];
    for (int ab = 0; ab < numPair; ++ab)
    {
        CONS_PAIR* consPair = constraint->conspairList[ab];
        int a = consPair->atomIindex;
        int b = consPair->atomJindex;
        VOP2(rab[ab], =, r[a], -, r[b]);
        nearestImage(&rab[ab].x, &rab[ab].y, &rab[ab].z);
    }
    double err = solveConstraintMatrix(location, numPair, numAtm, rab, constraint->conspairList, rMass, dt, lambda, v);
    int cnt = 1;
    while (err > tol)
    {
        err = solveConstraintMatrix(location, numPair, numAtm, rab, constraint->conspairList, rMass, dt, lambda, v);
        cnt++;
    }
    //printf("%d %e\n",cnt,err); 
    for (int j = 0; j < numAtm; ++j) // Store the new values
    {
        int index = resRangeStart + constraint->atomIDList[j];
        state->vx[index] = v[j].x;
        state->vy[index] = v[j].y;
        state->vz[index] = v[j].z;
    }
    //}
}

void paddingCons(double dt, STATE *state, NGLFCONSTRAINT_PARMS *p, GID_ORDER *gidOrder2, RESRANGE *resRange2)
{



    STATE *statechpad = &(p->charmmpot_parms->statechpad);
    if (getSize(0) == 1) statechpad = state;
    GID_ORDER *gidOrder = p->charmmpot_parms->gidOrder;
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    unsigned start = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        resRange2[i].start = start;
        start += residueList[i].size;
        resRange2[i].end = start;
    }

    if (state != statechpad)
    {
        statechpad->nlocal = state->nlocal;
        statechpad->nion = residueSet->molSize;

        unsigned nion = state->nion;

        for (int i = nion; i < residueSet->molSize; i++)
        {
            gidOrder2[i].id = INT_MAX; //this will be removed in next commit. see note below on why we set gid to this high number. 
            statechpad->label[i] = INT_MAX;
            ;
        }
        unsigned index = 0;
        unsigned s_id = 0;
        int sp_id = 0;
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
            char atmName[16];
            strcpy(atmName, name);
            strcat_(atmName, terminus);
            int len = strlen(atmName);
            for (int j = 0; j < resiConn->atomListSize; j++)
            {
                strcpy(atmName + len, resiConn->atomList[j]->atmName);
                s_id = gidOrder[index].id;
                if (index < nion && state->species[s_id] && (state->species[s_id]->name) &&(strcmp(atmName, state->species[s_id]->name) == 0))
                {
                    //if(statechpad->label[sp_id]!=state->label[s_id]) printf("pad gid=%llu  reg gid=%llu\n", statechpad->label[sp_id], state->label[s_id]);
                    gidOrder2[sp_id].id = gidOrder[index].id;
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
                }
                sp_id++;
            }
        }
    }
    //profile(CHARMM_COVALENT, END);
}

void copyState2Statechpad(STATE *state, STATE *statechpad, GID_ORDER *gidOrder2, SETLIST *residueSet)
{
    for (int i = 0; i < residueSet->molSize; i++)
    {
        int index = gidOrder2[i].id;
        if (index == INT_MAX) continue;
        if (index > state->nion)continue;
        statechpad->rx[i] = state->rx[index];
        statechpad->ry[i] = state->ry[index];
        statechpad->rz[i] = state->rz[index];
        statechpad->vx[i] = state->vx[index];
        statechpad->vy[i] = state->vy[index];
        statechpad->vz[i] = state->vz[index];
    }
}

void copyStatechpad2State(STATE *state, STATE *statechpad, GID_ORDER *gidOrder2, SETLIST *residueSet)
{
    for (int i = 0; i < residueSet->molSize; i++)
    {
        int index = gidOrder2[i].id;
        if (index == INT_MAX) continue;
        if (index > state->nion)continue;
        state->vx[index] = statechpad->vx[i];
        state->vy[index] = statechpad->vy[i];
        state->vz[index] = statechpad->vz[i];
    }
}

int pruneConsList(RESI_CONN *resiConn, GID_ORDER *gidOrder, int nlocal, int nTotal, int resRangeStart, CONSTRAINT **consList)
{
    int consListSize = 0;
    for (int i = 0; i < resiConn->consListSize; i++)
    {
        CONSTRAINT* constraint = resiConn->consList[i];
        //      int indexConsAtm=resRangeStart+constraint->atomID;
        //      if(indexConsAtm>nTotal)continue; //skip if atom doesn't exist
        //      if (gidOrder[indexConsAtm].id>= (int) nlocal)continue; //Not in local
        if (constraint->numPair == 0) continue;
        consList[consListSize++] = constraint;
    }
    return consListSize;
}

void velocityConstraintOld(double dt, STATE *state, NGLFCONSTRAINT_PARMS *p, GID_ORDER *gidOrder2, RESRANGE *resRange2, int location)
{
    STATE *statechpad = &p->charmmpot_parms->statechpad;
    if (getSize(0) == 1) statechpad = state;
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    if (state != statechpad) copyState2Statechpad(state, statechpad, gidOrder2, residueSet);
    LISTNODE* residueList = residueSet->list;
    THREE_SMATRIX sion = szero;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN *resiConn = residueList[i].resiConn;
        //CONSTRAINT *consList[resiConn->consListSize];
        //int consListSize  = pruneConsList(resiConn,gidOrder2,state->nlocal,residueSet->molSize,resRange2[i].start,consList);

        resMoveConsOld(resiConn->consListSize, resiConn->consList, dt, state, &sion, resRange2[i].start, location);
    }
    if (state != statechpad) copyStatechpad2State(state, statechpad, gidOrder2, residueSet);
}

void velocityConstraint(double dt, STATE *state, NGLFCONSTRAINT_PARMS *p, GID_ORDER *gidOrder2, RESRANGE *resRange2, int location)
{
    STATE *statechpad = &p->charmmpot_parms->statechpad;
    if (getSize(0) == 1) statechpad = state;
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    if (state != statechpad) copyState2Statechpad(state, statechpad, gidOrder2, residueSet);
    LISTNODE* residueList = residueSet->list;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN *resiConn = residueList[i].resiConn;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            CONSTRAINT* constraint = resiConn->consList[j];
            resMoveCons(constraint, dt, state, resRange2[i].start, location);
        }
    }
    if (state != statechpad) copyStatechpad2State(state, statechpad, gidOrder2, residueSet);
}

void velocityConstraintMolecule(double dt, STATE *state, NGLFCONSTRAINT_PARMS *p, GID_ORDER *gidOrder2, RESRANGE *resRange2, int location)
{
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN *resiConn = residueList[i].resiConn;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            CONSTRAINT* constraint = resiConn->consList[j];
            resMoveCons(constraint, dt, state, resRange2[i].start, location);
        }
    }
}

int countConstraints(NGLFCONSTRAINT_PARMS *p)
{
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    int nConstraints = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN *resiConn = residueList[i].resiConn;
        for (int i = 0; i < resiConn->consListSize; i++)
        {
            CONSTRAINT* constraint = resiConn->consList[i];
            //int numAtm = constraint->numAtom; // total number of atom
            int numPair = constraint->numPair;
            nConstraints += numPair;
        }
    }
    return nConstraints;
}

void nglfconstraint(DDC*ddc, SIMULATE*simulate, NGLFCONSTRAINT_PARMS*p)
{
    double dt = simulate->dt;

    SYSTEM* sys = simulate->system;
    STATE* state = sys->collection->state;
    double time = simulate->time;

    if (p->beta > 0)
    {
        THREE_SMATRIX pTensor = molecularPressure(sys, sys->energyInfo.virial, p->T);

        pTensor.xx -= p->P0;
        pTensor.yy -= p->P0;
        pTensor.zz -= p->P0;

        changeVolume(state, sys->box, &pTensor, p->beta, p->tauBarostat, dt);
    }

    for (int kk = 0; kk < state->nlocal; kk++) state->group[kk]->velocityUpdate(FRONT_TIMESTEP, kk, state->group[kk], state, time, 0.5 * dt);

    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    GID_ORDER gidOrder2[residueSet->molSize]; // Need to put on heap
    RESRANGE resRange2[residueSet->listSize]; // Need to put on heap
    paddingCons(dt, state, p, gidOrder2, resRange2);

    velocityConstraintOld(dt, state, p, gidOrder2, resRange2, FRONT_TIMESTEP);
    REF *v0 = vaf_v0();

    for (int kk = 0; kk < state->nlocal; kk++)
    {
        THREE_VECTOR delta;
        state->rx[kk] += delta.x = dt * state->vx[kk];
        state->ry[kk] += delta.y = dt * state->vy[kk];
        state->rz[kk] += delta.z = dt * state->vz[kk];
        if (v0 != NULL) VOP1(v0[kk].r, +=, delta); // for GPU code you could assume v0 == NULL;  can there VOP1 and delta are not needed. Shiv talk to me about this.  Jim 
    }

    for (int kk = 0; kk < state->nlocal; kk++) backInBox_fast(state->rx + kk, state->ry + kk, state->rz + kk);

    ddc->update = 0;
    time += dt; // positions, box (volume, h0, hinv), and forces at  t = n*dt + dt 
    simulate->time = sys->time = time;
    simulate->loop++;
    sys->loop = simulate->loop;

    for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update1(sys->group[kk], -1, state, time, 0.5 * dt);

    if (ddcenergy(ddc, sys, 0) != 0) return;

    for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update(sys->group[kk], BACK_TIMESTEP, state, time, 0.5 * dt);

    for (int kk = 0; kk < state->nlocal; kk++) state->group[kk]->velocityUpdate(BACK_TIMESTEP, kk, state->group[kk], state, time, 0.50 * dt);

    velocityConstraintOld(dt, state, p, gidOrder2, resRange2, BACK_TIMESTEP);

    kinetic_terms(sys, 1);

    for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update2(sys->group[kk], -1, state, time, 0.5 * dt);
    //eval_energyInfo(sys);
    for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update(sys->group[kk], FRONT_TIMESTEP, state, time, 0.5 * dt);

    simulate->time = sys->time = time;
    /*errorCheck(ddc->domain_id, simulate->loop, state, sys->energyInfo, p, datafile); */
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
