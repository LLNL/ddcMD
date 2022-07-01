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

THREE_VECTOR vector_sub(THREE_VECTOR a1, THREE_VECTOR a2)
{
    THREE_VECTOR a;
    a.x = a1.x - a2.x;
    a.y = a1.y - a2.y;
    a.z = a1.z - a2.z;
    return a;
}

THREE_VECTOR vector_add(THREE_VECTOR a1, THREE_VECTOR a2)
{
    THREE_VECTOR a;
    a.x = a1.x + a2.x;
    a.y = a1.y + a2.y;
    a.z = a1.z + a2.z;
    return a;
}

THREE_VECTOR vector_scalesub(double scale, THREE_VECTOR a1, THREE_VECTOR a2)
{
    THREE_VECTOR a;
    a.x = scale*a1.x - a2.x;
    a.y = scale*a1.y - a2.y;
    a.z = scale*a1.z - a2.z;
    return a;
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

void coor_vsite3fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;

    THREE_VECTOR r23;
    THREE_VECTOR r34;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);

    THREE_VECTOR r234=vector_sadd(a, r34, r23);
    double rsq = dot1(r234, r234);
    double rnom =sqrt(rsq);
    double b1 = b/rnom;

    state->rx[index1]= state->rx[index2] + b1 * r234.x;
    state->ry[index1]= state->ry[index2] + b1 * r234.y;
    state->rz[index1]= state->rz[index2] + b1 * r234.z;

    //Velocity
    THREE_VECTOR v23;
    THREE_VECTOR v34;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];

    v34.x = state->vx[index4] - state->vx[index3];
    v34.y = state->vy[index4] - state->vy[index3];
    v34.z = state->vz[index4] - state->vz[index3];

    THREE_VECTOR v234=vector_sadd(a, v34, v23);
    double rv = dot1(r234, v234);
    double c1 = rv/rsq;

    state->vx[index1]= state->vx[index2] + b1 * (v234.x - c1 * r234.x);
    state->vy[index1]= state->vy[index2] + b1 * (v234.y - c1 * r234.y);
    state->vz[index1]= state->vz[index2] + b1 * (v234.z - c1 * r234.z);
}

void coor_vsite3fad(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;

    THREE_VECTOR r23;
    THREE_VECTOR r34;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);

    double r23sq = dot1(r23, r23);
    double rnorm = sqrt(r23sq);
    double cr234 = dot1(r23, r34);
    double c1 = cr234/r23sq;
    THREE_VECTOR rp = vector_sadd(-1*c1, r23, r34);

    double a1 = a/rnorm;
    double rpsq = dot1(rp, rp);
    double b1 = b/ sqrt(rpsq);

    state->rx[index1]= state->rx[index2] + a1 * r23.x + b1 * rp.x;
    state->ry[index1]= state->ry[index2] + a1 * r23.y + b1 * rp.y;
    state->rz[index1]= state->rz[index2] + a1 * r23.z + b1 * rp.z;

    //velocity
    THREE_VECTOR v23;
    THREE_VECTOR v34;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];

    v34.x = state->vx[index4] - state->vx[index3];
    v34.y = state->vy[index4] - state->vy[index3];
    v34.z = state->vz[index4] - state->vz[index3];

    double cvr1 = dot1(v23, r34) + dot1(r23, v34);
    double crv2 = dot1(r23, v23);
    double a1v = cr234/r23sq;
    double a1r = (cvr1 - 2* a1v*crv2)/r23sq;

    THREE_VECTOR vtemp = vector_sadd(-1*a1r, r23, v34);
    THREE_VECTOR vp = vector_sadd(-1*a1v, v23, vtemp);
    double crvp = dot1(rp, vp);
    double a2r = crv2/r23sq;
    double b2rv = crvp/rpsq;

    THREE_VECTOR va2r = vector_sadd(-1*a2r, r23, v23);
    THREE_VECTOR vb2rv = vector_sadd(-1*b2rv, rp, vp);

    state->vx[index1]= state->vx[index2] + a1 * va2r.x + b1 * vb2rv.x;
    state->vy[index1]= state->vy[index2] + a1 * va2r.y + b1 * vb2rv.y;
    state->vz[index1]= state->vz[index2] + a1 * va2r.z + b1 * vb2rv.z;
}

void coor_vsite4fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5){
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r34;
    THREE_VECTOR r35;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);
    bioVec(state, index5, index3, &r35);

    THREE_VECTOR rcom1 = vector_sadd(a, r34, r23);
    THREE_VECTOR rcom = vector_sadd(b, r35, rcom1);

    double rcomsq = dot1(rcom, rcom);
    double c1 = c /sqrt(rcomsq);

    state->rx[index1]= state->rx[index2] + c1 * rcom.x;
    state->ry[index1]= state->ry[index2] + c1 * rcom.y;
    state->rz[index1]= state->rz[index2] + c1 * rcom.z;

    //velocity
    THREE_VECTOR v23;
    THREE_VECTOR v34;
    THREE_VECTOR v35;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];

    v34.x = state->vx[index4] - state->vx[index3];
    v34.y = state->vy[index4] - state->vy[index3];
    v34.z = state->vz[index4] - state->vz[index3];

    v35.x = state->vx[index5] - state->vx[index3];
    v35.y = state->vy[index5] - state->vy[index3];
    v35.z = state->vz[index5] - state->vz[index3];

    THREE_VECTOR vcom1 = vector_sadd(a, v34, v23);
    THREE_VECTOR vcom = vector_sadd(b, v35, vcom1);

    double cvr = dot1(vcom, rcom);
    double c1v = cvr/rcomsq;
    THREE_VECTOR vrcom = vector_sadd(-1*c1v, rcom, vcom);

    state->vx[index1]= state->vx[index2] + c1 * vrcom.x;
    state->vy[index1]= state->vy[index2] + c1 * vrcom.y;
    state->vz[index1]= state->vz[index2] + c1 * vrcom.z;
}

void coor_vsite4fdn(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5){
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r24;
    THREE_VECTOR r25;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);
    bioVec(state, index5, index2, &r25);

    THREE_VECTOR rja = vector_scalesub(a, r24, r23);
    THREE_VECTOR rjb = vector_scalesub(b, r25, r23);

    THREE_VECTOR rjcom = cross(&rja, &rjb);
    double rjcomsq = dot1(rjcom, rjcom);
    double c1 = c / sqrt(rjcomsq);

    state->rx[index1]= state->rx[index2] + c1 * rjcom.x;
    state->ry[index1]= state->ry[index2] + c1 * rjcom.y;
    state->rz[index1]= state->rz[index2] + c1 * rjcom.z;

    //velocity
    THREE_VECTOR v23;
    THREE_VECTOR v24;
    THREE_VECTOR v25;

    v23.x = state->vx[index3] - state->vx[index2];
    v23.y = state->vy[index3] - state->vy[index2];
    v23.z = state->vz[index3] - state->vz[index2];

    v24.x = state->vx[index4] - state->vx[index2];
    v24.y = state->vy[index4] - state->vy[index2];
    v24.z = state->vz[index4] - state->vz[index2];

    v25.x = state->vx[index5] - state->vx[index2];
    v25.y = state->vy[index5] - state->vy[index2];
    v25.z = state->vz[index5] - state->vz[index2];

    THREE_VECTOR vja;
    THREE_VECTOR vjb;
    vja.x = a * v24.x - v23.x;
    vja.y = a * v24.y - v23.y;
    vja.z = a * v24.z - v23.z;
    vjb.x = b * v25.x - v23.x;
    vjb.y = b * v25.y - v23.y;
    vjb.z = b * v25.z - v23.z;

    THREE_VECTOR vcom1 = cross(&vja, &rjb);
    THREE_VECTOR vcom2 = cross(&rja, &vjb);
    THREE_VECTOR vcom;
    vcom.x = vcom1.x + vcom2.x;
    vcom.y = vcom1.y + vcom2.y;
    vcom.z = vcom1.z + vcom2.z;

    double cvrcom = dot1(rjcom, vcom);
    double c1v = cvrcom/rjcomsq;
    THREE_VECTOR vrcom = vector_sadd(-1*c1v, rjcom, vcom);

    state->vx[index1]= state->vx[index2] + c1 * vrcom.x;
    state->vy[index1]= state->vy[index2] + c1 * vrcom.y;
    state->vz[index1]= state->vz[index2] + c1 * vrcom.z;
}

void coor_vsiten(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5){
    double a=vsiteConn->a;
    THREE_VECTOR r23;
    THREE_VECTOR r24;
    THREE_VECTOR r25;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);
    bioVec(state, index5, index2, &r25);

    state->rx[index1]= state->rx[index2] + a * (r23.x + r24.x + r25.x);
    state->ry[index1]= state->ry[index2] + a * (r23.y + r24.y + r25.y);
    state->rz[index1]= state->rz[index2] + a * (r23.z + r24.z + r25.z);

    //velocity
    state->vx[index1]= a * (state->vx[index2] + state->vx[index3] + state->vx[index4] + state->vx[index5]);
    state->vy[index1]= a * (state->vy[index2] + state->vy[index3] + state->vy[index4] + state->vy[index5]);
    state->vz[index1]= a * (state->vz[index2] + state->vz[index3] + state->vz[index4] + state->vz[index5]);
}


void vsite_resdiue_construction(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms)
{
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    RESI_CONN* resiConn = findResiConnNew(charmmParms, name);

    //for (int i = 0; i < resiConn->vsiteListSize; ++i) {
    //    VSITE_CONN *vsiteConn=resiConn->vsiteList[i];
    VSITE_CONN **vsiteList = resiConn->vsiteList;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE vsiteRange = resiConn->atmRanges[i]->vsiteRange;
        if (vsiteRange.start == -1)
        {
            continue;
        }
        for (int vsiteIndex = vsiteRange.start; vsiteIndex < vsiteRange.end; vsiteIndex++) {
            VSITE_CONN *vsiteConn = vsiteList[vsiteIndex];
            int index1 = resRange->start + vsiteConn->atom1;
            int index2 = resRange->start + vsiteConn->atom2;
            int index3 = resRange->start + vsiteConn->atom3;
            int index4 = resRange->start + vsiteConn->atom4;
            int index5 = resRange->start + vsiteConn->atom5;

            switch (vsiteConn->vtype) {
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
                case VSITE3FD:
                    coor_vsite3fd(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE3FAD:
                    coor_vsite3fad(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE3OUT:
                    coor_vsite3out(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE4FD:
                    coor_vsite4fd(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                case VSITE4FDN:
                    coor_vsite4fdn(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                case VSITEN:
                    coor_vsiten(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                default:
                    printf("Virtual site type %d is not support", vsiteConn->vtype);
            }
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

void force_vsite3fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    THREE_VECTOR r23;
    THREE_VECTOR r34;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);

    THREE_VECTOR r234=vector_sadd(a, r34, r23);
    double rsq = dot1(r234, r234);
    double rnom =sqrt(rsq);
    double b1 = b/rnom;
    double b2 = dot1(r234, f1)/rsq;

    THREE_VECTOR rcom = vector_sadd(-1*b2, r234, f1);
    VSCALE(rcom, b1);

    double a1 = 1-a;
    THREE_VECTOR f2 = vector_sadd(-1, rcom, f1);
    accufVec(state, index2, f2);
    THREE_VECTOR f3 = rcom;
    VSCALE(f3, a1);
    accufVec(state, index3, f3);
    THREE_VECTOR f4 = rcom;
    VSCALE(f4, a);
    accufVec(state, index4, f4);

}

void force_vsite3fad(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    THREE_VECTOR r23;
    THREE_VECTOR r34;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);

    double r23sq = dot1(r23, r23);
    double r23norm = sqrt(r23sq);
    double c1 = dot1(r23, r34)/r23sq;

    THREE_VECTOR rp = vector_sadd(-1*c1, r23, r34);
    double rpsq = dot1(rp, rp);
    double rpnorm = sqrt(rpsq);
    double a1 = a/r23norm;
    double b1 = b/rpnorm;

    double cfr1 = dot1(r23, f1)/r23sq;
    THREE_VECTOR fpij = r23;
    VSCALE(fpij, cfr1);

    double cfr2 = dot1(rp, f1)/rpsq;
    THREE_VECTOR fppp = rp;
    VSCALE(fppp, cfr2);

    THREE_VECTOR f2 = vector_sub(f1, fpij);
    VSCALE(f2, a1);

    THREE_VECTOR f3 = vector_sub(f2, fppp);
    VSCALE(f3, b1);

    THREE_VECTOR f4 = rp;
    VSCALE(f4, b1 * cfr1);

    THREE_VECTOR f2com1 = vector_sub(f1, f2);
    THREE_VECTOR f2com2 = vector_sadd(c1, f3, f2com1);
    THREE_VECTOR f2com = vector_add(f2com2, f4);
    accufVec(state, index2, f2com);

    double c3 = c1 +1;
    THREE_VECTOR f3com1 = vector_sadd(-1* c3, f3, f2);
    THREE_VECTOR f3com = vector_sub(f3com1, f4);
    accufVec(state, index3, f3com);

    accufVec(state, index4, f3);
}

void force_vsite4fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r34;
    THREE_VECTOR r35;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index3, &r34);
    bioVec(state, index5, index3, &r35);

    THREE_VECTOR rcom1 = vector_sadd(a, r34, r23);
    THREE_VECTOR rcom = vector_sadd(b, r35, rcom1);

    double rcomsq = dot1(rcom, rcom);
    double c1 = c /sqrt(rcomsq);

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    double crf = dot1(rcom, f1)/rcomsq;
    THREE_VECTOR fcom = (-1*crf, rcom, f1);
    VSCALE(fcom, c1);

    double a1 = 1 - a - b ;
    THREE_VECTOR f2 = vector_sub(f1, fcom);
    accufVec(state, index2, f2);

    THREE_VECTOR f3 = fcom;
    VSCALE(f3, a1);
    accufVec(state, index3, f3);

    THREE_VECTOR f4 = fcom;
    VSCALE(f4, a);
    accufVec(state, index4, f4);

    THREE_VECTOR f5 = fcom;
    VSCALE(f5, b);
    accufVec(state, index5, f5);
}

void force_vsite4fdn(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5)
{
    double a=vsiteConn->a;
    double b=vsiteConn->b;
    double c=vsiteConn->c;

    THREE_VECTOR r23;
    THREE_VECTOR r24;
    THREE_VECTOR r25;

    bioVec(state, index3, index2, &r23);
    bioVec(state, index4, index2, &r24);
    bioVec(state, index5, index2, &r25);

    THREE_VECTOR rja = vector_scalesub(a, r24, r23);
    THREE_VECTOR rjb = vector_scalesub(b, r25, r23);
    THREE_VECTOR rab = vector_sub(rjb, rja);

    THREE_VECTOR rjcom = cross(&rja, &rjb);
    double rjcomsq = dot1(rjcom, rjcom);
    double rjnorm = sqrt(rjcomsq);
    double c1 = c / rjnorm;

    THREE_VECTOR f1;
    scalefVec(state, index1, 1.0, &f1);

    THREE_VECTOR fcom=f1;
    VSCALE(fcom, c/rjnorm);

    THREE_VECTOR rt = cross(&rjcom, &rab);
    VSCALE(rt, 1/rjcomsq);

    THREE_VECTOR f3;
    f3.x = (-rjcom.x * rt.x) * fcom.x + (rab.z - rjcom.y * rt.x) * fcom.y + (-rab.y - rjcom.z * rt.x) * fcom.z;
    f3.y = (-rab.z - rjcom.x * rt.y) * fcom.x + (-rjcom.y * rt.y) * fcom.y + (rab.x - rjcom.z * rt.y) * fcom.z;
    f3.z = (rab.y - rjcom.x * rt.z) * fcom.x + (-rab.x -rjcom.y * rt.z) * fcom.y + (-rjcom.z * rt.z) * fcom.z;

    rt = cross(&rjb, &rjcom);
    VSCALE(rt, a/rjcomsq);
    THREE_VECTOR f4;
    f4.x = (-rjcom.x * rt.x) * fcom.x + (-a * rjb.z -rjcom.y * rt.x) * fcom.y + (a * rjb.y - rjcom.z * rt.z) * fcom.z;
    f4.y = (a * rjb.z - rjcom.x * rt.y) * fcom.x + (-rjcom.y * rt.y) * fcom.y + (-a * rjb.x - rjcom.z * rt.y) * fcom.z;
    f4.z = (-a * rjb.y - rjcom.x * rt.z) * fcom.x + (a * rjb.x - rjcom.x * rt.z) * fcom.y + (-rjcom.z * rt.z) * fcom.z;

    rt = cross(&rjcom, &rja);
    VSCALE(rt, b/rjcomsq);
    THREE_VECTOR f5;
    f5.x = (-rjcom.x * rt.x) * fcom.x + (b * rja.z - rjcom.y * rt.x) * fcom.y + (-b * rja.y - rjcom.z * rt.x) * fcom.z;
    f5.y = (-b * rja.z - rjcom.x * rt.y) * fcom.x + (-rjcom.y * rt.y) * fcom.y + (b * rja.x - rjcom.z * rt.y) * fcom.z;
    f5.y = (b * rja.y - rjcom.x * rt.z) * fcom.x + (-b * rja.x - rjcom.y * rt.z) * fcom.y + (-rjcom.z * rt.z) * fcom.z;

    THREE_VECTOR f2com1 = vector_sub(f1, f3);
    THREE_VECTOR f2com2 = vector_sub(f2com1, f4);
    THREE_VECTOR f2com  = vector_sub(f2com2, f5);
    accufVec(state, index2, f2com);

    accufVec(state, index3, f3);
    accufVec(state, index4, f4);
    accufVec(state, index5, f5);
}

void force_vsiten(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4, int index5)
{
    double a=vsiteConn->a;

    THREE_VECTOR f1;
    f1.x = a* state->fx[index1];
    f1.y = a* state->fy[index1];
    f1.z = a* state->fz[index1];

    accufVec(state, index2, f1);
    accufVec(state, index3, f1);
    accufVec(state, index4, f1);
    accufVec(state, index5, f1);

}

void vsite_residue_force(STATE* state, GID_ORDER* gidOrder, unsigned nlocal, unsigned nTotal, char * name, RESRANGE* resRange, CHARMMPOT_PARMS *parms) {
    CHARMM_PARMS* charmmParms = parms->charmmParms;
    RESI_CONN* resiConn = findResiConnNew(charmmParms, name);

    //for (int i = 0; i < resiConn->vsiteListSize; ++i) {
    //    VSITE_CONN *vsiteConn=resiConn->vsiteList[i];
    VSITE_CONN **vsiteList = resiConn->vsiteList;
    for (int i = 0; i < resiConn->atomListSize; ++i)
    {
        RANGE vsiteRange = resiConn->atmRanges[i]->vsiteRange;
        if (vsiteRange.start == -1)
        {
            continue;
        }
        for (int vsiteIndex = vsiteRange.start; vsiteIndex < vsiteRange.end; vsiteIndex++) {
            VSITE_CONN* vsiteConn = vsiteList[vsiteIndex];
            int index1 = resRange->start + vsiteConn->atom1;
            int index2 = resRange->start + vsiteConn->atom2;
            int index3 = resRange->start + vsiteConn->atom3;
            int index4 = resRange->start + vsiteConn->atom4;
            int index5 = resRange->start + vsiteConn->atom5;

            switch (vsiteConn->vtype) {
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
                case VSITE3FD:
                    force_vsite3fd(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE3FAD:
                    force_vsite3fad(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE3OUT:
                    force_vsite3out(state, vsiteConn, index1, index2, index3, index4);
                    break;
                case VSITE4FD:
                    force_vsite4fd(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                case VSITE4FDN:
                    force_vsite4fdn(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                case VSITEN:
                    force_vsiten(state, vsiteConn, index1, index2, index3, index4, index5);
                    break;
                default:
                    printf("Virtual site type %d is not support", vsiteConn->vtype);
            }
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