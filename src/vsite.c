#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "collection.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"

void vsite_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in)
{
}
void vsite_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt)
{
    /*
    double *vx = state->vx;
    double *vy = state->vy;
    double *vz = state->vz;
    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;
    SPECIES **species = state->species;
    double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;

    double a = dt/mass;
    vx[k] += a*fx[k] ;
    vy[k] += a*fy[k] ;
    vz[k] += a*fz[k] ;
     */
}

void vsite_parms(GROUP *gp)
{
    gp->itype = FREE;
    gp->parm = NULL;
    gp->write_dynamics = NULL;
    gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))vsite_velocityUpdate;
    gp->Update= (void (*)(GROUP *, int, void *, double, double))vsite_Update;
}

void coor_vsite1(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2)
{
    state->rx[index1]= state->rx[index2];
    state->ry[index1]= state->ry[index2];
    state->rz[index1]= state->rz[index2];

    state->vx[index1]= state->vx[index2];
    state->vy[index1]= state->vy[index2];
    state->vz[index1]= state->vz[index2];
}

void coor_vsite2(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3)
{
    double a=vsiteConn->a;
    double b= 1- a;

    THREE_VECTOR r23;
    bioVec(state, index3, index2, &r23);

    state->rx[index1]= state->rx[index2] + a * r23.x;
    state->ry[index1]= state->ry[index2] + a * r23.y;
    state->rz[index1]= state->rz[index2] + a * r23.z;

    // Velocity
    state->vx[index1]= b * state->vx[index2] + a * state->vx[index3];
    state->vy[index1]= b * state->vy[index2] + a * state->vy[index3];
    state->vz[index1]= b * state->vz[index2] + a * state->vz[index3];
}

void coor_vsite2fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3)
{
    double a=vsiteConn->a;

    THREE_VECTOR r23;
    bioVec(state, index3, index2, &r23);
    double invRsq=1.0/dot1(r23, r23);
    double b= a*sqrt(invRsq);

    state->rx[index1]= state->rx[index2] + b * r23.x;
    state->ry[index1]= state->ry[index2] + b * r23.y;
    state->rz[index1]= state->rz[index2] + b * r23.z;

    // Velocity
    THREE_VECTOR v23;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];
    double v23r23 = dot1(v23, r23);
    double ceof = v23r23 * invRsq;

    state->vx[index1]= state->vx[index2] + b * (v23.x - ceof * r23.x);
    state->vy[index1]= state->vy[index2] + b * (v23.y - ceof * r23.y);
    state->vz[index1]= state->vz[index2] + b * (v23.z - ceof * r23.z);

}

void coor_vsite3(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c= 1 - a - b;

    THREE_VECTOR r23;
    THREE_VECTOR r24;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);

    state->rx[index1]= state->rx[index2] + a * r23.x + b * r24.x;
    state->ry[index1]= state->ry[index2] + a * r23.y + b * r24.y;
    state->rz[index1]= state->rz[index2] + a * r23.z + b * r24.z;

    // Velocity
    state->vx[index1]= c * state->vx[index2] + a * state->vx[index3] + b * state->vx[index4];
    state->vy[index1]= c * state->vy[index2] + a * state->vy[index3] + b * state->vy[index4];
    state->vz[index1]= c * state->vz[index2] + a * state->vz[index3] + b * state->vz[index4];

}

void coor_vsite3out(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r24;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);

    THREE_VECTOR rx = cross(&r23, &r24);

    state->rx[index1]= state->rx[index2] + a * r23.x + b * r24.x + c * rx.x;
    state->ry[index1]= state->ry[index2] + a * r23.y + b * r24.y + c * rx.y;
    state->rz[index1]= state->rz[index2] + a * r23.z + b * r24.z + c * rx.z;

    //Velocity
    THREE_VECTOR v23;
    THREE_VECTOR v24;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];

    v24.x = state->vx[index4] - state->vx[index2];
    v24.y = state->vy[index4] - state->vy[index2];
    v24.z = state->vz[index4] - state->vz[index2];

    THREE_VECTOR v23r24 = cross(&v23, &r24);
    THREE_VECTOR r23v24 = cross(&r23, &v24);
    state->vx[index1]= state->vx[index2] + a * v23.x + b * v24.x + c * (v23r24.x + r23v24.x);
    state->vy[index1]= state->vy[index2] + a * v23.y + b * v24.y + c * (v23r24.y + r23v24.y);
    state->vz[index1]= state->vz[index2] + a * v23.z + b * v24.z + c * (v23r24.z + r23v24.z);
}

void vsite_resdiue_construction(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms)
{
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    RESI_CONN* resiConn = findResiConnNew(charmmParms, name);

    for (int i = 0; i < resiConn->vsiteListSize; ++i) {
        VSITE_CONN *vsiteConn=resiConn->vsiteList[i];
        int index1 = resRange->start + vsiteConn->atom1;
        int index2 = resRange->start + vsiteConn->atom2;
        int index3 = resRange->start + vsiteConn->atom3;
        int index4 = resRange->start + vsiteConn->atom4;

        switch(vsiteConn->vtype)
        {
            case VSITE1:
                coor_vsite1(state, vsiteConn, index1, index2);
                break;
            case VSITE2:
                coor_vsite2(state, vsiteConn, index1, index2, index3);
                break;
            case VSITE2FD:
                coor_vsite2fd(state, vsiteConn, index1, index2, index3);
                break;
            case VSITE3:
                coor_vsite3(state, vsiteConn, index1, index2, index3, index4);
                break;
            case VSITE3OUT:
                coor_vsite3out(state, vsiteConn, index1, index2, index3, index4);
                break;
        }

    }

}

void vsite_construction(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    SETLIST *residueSet = &parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    STATE *statechpad = &(parms->statechpad);

    for (int i = 0; i < residueSet->listSize; i++)
    {
        vsite_resdiue_construction(statechpad, parms->gidOrder2, sys->nlocal, residueSet->molSize, residueList[i].name, &(parms->resRange2[i]), parms);
    }
}

void scalefVec(STATE* state, int index, double coef, THREE_VECTOR* vec)
{
    vec->x = coef*state->fx[index];
    vec->y = coef*state->fy[index];
    vec->z = coef*state->fz[index];
}

void accufVec(STATE* state, int index, THREE_VECTOR vec)
{
    state->fx[index]+=vec.x;
    state->fy[index]+=vec.y;
    state->fz[index]+=vec.z;
}

void force_vsite1(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2)
{
    state->fx[index2] += state->fx[index1];
    state->fy[index2] += state->fy[index1];
    state->fz[index2] += state->fz[index1];
}

void force_vsite2(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3)
{
    THREE_VECTOR f2;
    THREE_VECTOR f3;

    double a=vsiteConn->a;
    scalefVec(state, index1, 1-a, &f2);
    scalefVec(state, index1, a, &f3);

    accufVec(state, index2, f2);
    accufVec(state, index3, f3);

}

void force_vsite2fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3)
{
    double a=vsiteConn->a;

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    THREE_VECTOR r23;
    bioVec(state, index3, index2, &r23);

    double invRsq=1.0/dot1(r23, r23);
    double b= a*sqrt(invRsq);
    double ceof = -1 * dot1(r23, f1) * invRsq;

    THREE_VECTOR f3 = vector_sadd(ceof, r23, f1);
    VSCALE(f3, b);

    state->fx[index2]+=f1.x - f3.x;
    state->fy[index2]+=f1.y - f3.y;
    state->fz[index2]+=f1.z - f3.z;

    accufVec(state, index3, f3);
}

void force_vsite3(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    THREE_VECTOR f2;
    THREE_VECTOR f3;
    THREE_VECTOR f4;

    double a=vsiteConn->a;
    double b=vsiteConn->b;

    scalefVec(state, index1, 1-a-b, &f2);
    scalefVec(state, index1, a, &f3);
    scalefVec(state, index1, b, &f4);

    accufVec(state, index2, f2);
    accufVec(state, index3, f3);
    accufVec(state, index4, f4);
}

void force_vsite3out(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r24;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    THREE_VECTOR cf1;
    scalefVec(state, index1, c, &cf1);

    THREE_VECTOR f3;
    THREE_VECTOR f4;

    f3.x = a * f1.x - r24.z * cf1.y + r24.y * cf1.z;
    f3.y = r24.z * cf1.x + a * f1.y - r24.x * cf1.z;
    f3.z = - r24.y * cf1.x + r24.x * cf1.y + a * f1.z;

    f4.x = b * f1.x + r23.z * cf1.y - r23.y * cf1.z;
    f4.y = -r23.z * cf1.x + b * f1.y + r23.x * cf1.z;
    f4.z = r23.y * cf1.x -r23.x * cf1.y + b * f1.z;

    state->fx[index2]+=f1.x - f3.x - f4.x;
    state->fy[index2]+=f1.y - f3.y - f4.y;
    state->fz[index2]+=f1.z - f3.z - f4.z;

    accufVec(state, index3, f3);
    accufVec(state, index4, f4);

}

void vsite_residue_force(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms) {
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    RESI_CONN* resiConn = findResiConnNew(charmmParms, name);

    for (int i = 0; i < resiConn->vsiteListSize; ++i) {
        VSITE_CONN *vsiteConn=resiConn->vsiteList[i];
        int index1 = resRange->start + vsiteConn->atom1;
        int index2 = resRange->start + vsiteConn->atom2;
        int index3 = resRange->start + vsiteConn->atom3;
        int index4 = resRange->start + vsiteConn->atom4;

        switch(vsiteConn->vtype)
        {
            case VSITE1:
                force_vsite1(state, vsiteConn, index1, index2);
                break;
            case VSITE2:
                force_vsite2(state, vsiteConn, index1, index2, index3);
                break;
            case VSITE2FD:
                force_vsite2fd(state, vsiteConn, index1, index2, index3);
                break;
            case VSITE3:
                force_vsite3(state, vsiteConn, index1, index2, index3, index4);
                break;
            case VSITE3OUT:
                force_vsite3out(state, vsiteConn, index1, index2, index3, index4);
                break;
        }

    }

}

void vsite_force(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    SETLIST *residueSet = &parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    STATE *statechpad = &(parms->statechpad);

    for (int i = 0; i < residueSet->listSize; i++)
    {
        vsite_residue_force(statechpad, parms->gidOrder2, sys->nlocal, residueSet->molSize, residueList[i].name, &(parms->resRange2[i]), parms);
    }
}