#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "collection.h"

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