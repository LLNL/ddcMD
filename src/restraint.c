#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "object.h"
#include "ddcMalloc.h"
#include "restraint.h"
#include "system.h"
#include "gid.h"
#include "error.h"
#include "simulate.h"
#include "units.h"
#include "preduce.h"
#include "utilities.h"
#include "preduce.h"
#include "accelerator.h"

int compareRestraint(const void **v1, const void **v2)
{
    const RESTRAINTPARMS *u1 = *(RESTRAINTPARMS **) v1;
    const RESTRAINTPARMS *u2 = *(RESTRAINTPARMS **) v2;
    if (u1->gid < u2->gid)
        return -1;
    else if (u1->gid > u2->gid)
        return 1;
    else
        return 0;
}

RESTRAINTPARMS *restraintparms_init(void *parent, char *name)
{
    RESTRAINTPARMS *restraintparms;
    restraintparms = (RESTRAINTPARMS *) object_initialize(name, "RESTRAINTPARMS", sizeof (RESTRAINTPARMS));
    restraintparms->parent = parent;
    object_get((OBJECT *) restraintparms, "gid", &(restraintparms->gid), U64, 1, "0");
    object_get((OBJECT *) restraintparms, "atomI", &(restraintparms->atomI), INT, 1, "0");
    object_get((OBJECT *) restraintparms, "func", &(restraintparms->func), INT, 1, "1");
    object_get((OBJECT *) restraintparms, "fcx", &(restraintparms->fcx), INT, 1, "0");
    object_get((OBJECT *) restraintparms, "fcy", &(restraintparms->fcy), INT, 1, "0");
    object_get((OBJECT *) restraintparms, "fcz", &(restraintparms->fcz), INT, 1, "0");
    object_get((OBJECT *) restraintparms, "x0", &(restraintparms->x0), DOUBLE, 1, "0");
    object_get((OBJECT *) restraintparms, "y0", &(restraintparms->y0), DOUBLE, 1, "0");
    object_get((OBJECT *) restraintparms, "z0", &(restraintparms->z0), DOUBLE, 1, "0");
    object_get((OBJECT *) restraintparms, "kb", &(restraintparms->kb), WITH_UNITS, 1, "0.0", "kJ*mol^-1*nm^-2", NULL);

    return restraintparms;

}

RCUT_TYPE *restraintCutoff(SYSTEM* sys, void*parms, int *n)
{
    static RCUT_TYPE rcut[1];
    *n = 0;
    rcut[0].value = 0.0;
    rcut[0].mode = RCUT_ALL;

    return rcut;
}

RESTRAINTLIST *restraint_init(void *parent, char *name)
{

    RESTRAINTLIST * restraint = (RESTRAINTLIST *) object_initialize(name, "RESTRAINTLIST", sizeof (RESTRAINTLIST));
    restraint->parent = parent;
    object_get((OBJECT *) restraint, "origin", &(restraint->origin), INT, 1, "0");
    char** restraintNames;
    restraint->nRestraint = object_getv((OBJECT *) restraint, "restraintList", (void *) &restraintNames, STRING, IGNORE_IF_NOT_FOUND); // bond is not required.
    restraint->restraintList = (RESTRAINTPARMS**) ddcMalloc(restraint->nRestraint * sizeof (RESTRAINTPARMS*));
    for (int i = 0; i < restraint->nRestraint; i++)
    {
        //printf("%s\n", bondNames[i]);
        restraint->restraintList[i] = restraintparms_init(restraint, restraintNames[i]);
    }

    // Currently the x0, y0, z0 read in from restraint.data are relatvie position
    // need to recale by the box size
    SYSTEM *sys = (SYSTEM*) ((POTENTIAL*) parent)->parent;
    double xBox = sys->box->h0.xx;
    double yBox = sys->box->h0.yy;
    double zBox = sys->box->h0.zz;
    //THREE_MATRIX h0;
    //box_get(NULL, HO_PTR, (void *) &h0); // NULL to trick to get the current box
    double length_convert = units_convert(1.0, NULL, "Angstrom");
    char message[80];
    sprintf(message, "Current Box Size x =%10.3f y =%10.3f z =%10.3f",
            xBox*length_convert, yBox*length_convert, zBox * length_convert);
    timestamp(message);

    //sort the gid in the restraint list  
    /*
       RESTRAINTPARMS* restraintArray=*(restraint->restraintList);
       for(int i = 0; i < restraint->nRestraint; i++){
       RESTRAINTPARMS *restraintparms=restraint->restraintList[i];
       printf("Before gid=%llu atomI=%d fcz=%d fcx=%d z0=%f\n", 
       restraintparms->gid, restraintparms->atomI, restraintparms->fcz, restraintparms->fcx, restraintparms->z0);
       }
     */

#if __APPLE__
    qsort(restraint->restraintList, restraint->nRestraint, sizeof (RESTRAINTPARMS*),(int(*)(const void*,const void*))compareRestraint);
#else
    qsort(restraint->restraintList, restraint->nRestraint, sizeof (RESTRAINTPARMS*), (__compar_fn_t) compareRestraint);
#endif
    /*
       for(int i = 0; i < restraint->nRestraint; i++){
       RESTRAINTPARMS *restraintparms=restraint->restraintList[i];
       printf("After gid=%llu atomI=%d fcz=%d fcx=%d z0=%f\n", 
       restraintparms->gid, restraintparms->atomI, restraintparms->fcz, restraintparms->fcx, restraintparms->z0);
       }    
     */

    restraint->restGidList = (gid_type*) ddcMalloc(restraint->nRestraint * sizeof (gid_type));
    ;

    // Get a temporary list of restraint gid list
    for (int i = 0; i < restraint->nRestraint; i++)
    {
        restraint->restGidList[i] = restraint->restraintList[i]->gid;
    }

    return restraint;
}

int binarySearchGidListJim(gid_type queryGid, gid_type* gidList, int first, int last)
{

    if (first <= last)
    {
        int mid = (first + last) / 2;
        if (queryGid == gidList[mid]) return mid;
        else if (queryGid < gidList[mid])
            binarySearchGidListJim(queryGid, gidList, first, mid - 1);
        else
            binarySearchGidListJim(queryGid, gidList, mid + 1, last);
    }
    return -1;
}

int binarySearchGidList(gid_type queryGid, gid_type* gidList, int first, int last)
{

    if (first <= last)
    {
        int mid = (first + last) / 2;
        //printf("gid=%llu gidList[%d]=%llu f=%d l=%d \n", queryGid, mid, gidList[mid], first, last);
        if (queryGid == gidList[mid])
        {
            return mid;
        }
        else if (queryGid < gidList[mid])
        {
            return binarySearchGidList(queryGid, gidList, first, mid - 1);
        }
        else
        {
            return binarySearchGidList(queryGid, gidList, mid + 1, last);
        }
    }
    else
    {
        return -1;
    }
}

void assignRestraintMap(SYSTEM *sys, RESTRAINTLIST *restraint)
{
    int nRestraint = restraint->nRestraint;

    for (int i = 0; i < nRestraint; i++)
    {
        restraint->restraintMap[i] = -1;
    }

    STATE * state = sys->collection->state;

    for (unsigned int ii = 0; ii < sys->nlocal; ii++)
    { // nlocal in state has not been initialized yet
        //printf("gid=%llu\n", state->label[ii]);

        int restIndex = binarySearchGidList(state->label[ii], restraint->restGidList, 0, nRestraint - 1);
        if (restIndex >= 0)
        {
            restraint->restraintMap[restIndex] = ii;
            //           printf("%s ii=%d restIndex=%d\n", state->species[ii]->name, ii, restIndex);
        }
    }

    //  for(int i=0; i<nRestraint; i++){
    //      printf("r=%d     ii=%d\n", i, restraint->restraintMap[i]);
    //  }    

    //for(int i=0; i<nRestraint; i++){
    //    if(restraint->restraintMap[i]==-1){
    //        error_action("restraint_parms error: constraintMap is not fully assigned", ERROR_IN("restraint_parms", ABORT));
    //    }
    //}   


}

RESTRAINTLIST *restraint_parms(POTENTIAL *potential)
{
    timestamp("Start Restraint Setup");

    char *parmfile;
    object_get((OBJECT *) potential, "parmfile", &parmfile, STRING, 1, "restraint.data");

    object_compilefile(parmfile);

    RESTRAINTLIST *restraint = restraint_init(potential, "restraint");

    // Create map constraint to the state

    restraint->restraintMap = (unsigned *) ddcMalloc(restraint->nRestraint * sizeof (unsigned*));

    SYSTEM *sys = (SYSTEM*) potential->parent;

    assignRestraintMap(sys, restraint);

    ACCELERATOR *accelerator = NULL;
    accelerator = accelerator_getAccelerator(accelerator);
    if ((accelerator != NULL) && accelerator->itype == GPU_CUDA)
    {
        printf("using gpu martini parms\n");
        restraintGPU_parms(restraint);
        potential->eval_potential = (void (*)(void *, void *, void *))restraintGPU;
    }

    timestamp("END   Restraint Setup");
    return restraint;

}

void scaleRestraintByBoxChange(RESTRAINTLIST *parms)
{

    THREE_MATRIX hfac;
    THREE_VECTOR rold, r;
    box_get(NULL, HFAC, (void *) &hfac); // NULL to trick to get the current box
    //printf("Outside box change\n");
    if (!matrix_equal(hfac, I_3x3))
    {
        //printf("Inside box change\n");
        for (int i = 0; i < parms->nRestraint; i++)
        {
            RESTRAINTPARMS *restraintparms = parms->restraintList[i];

            rold.x = restraintparms->x0;
            rold.y = restraintparms->y0;
            rold.z = restraintparms->z0;
            r = matrix_vector(hfac, rold);
            restraintparms->x0 = r.x;
            restraintparms->y0 = r.y;
            restraintparms->z0 = r.z;

        }
    }
}

void restraint(SYSTEM*sys, RESTRAINTLIST *parms, ETYPE *e)
{
    STATE * state = sys->collection->state;
    SIMULATE* simulation = (SIMULATE*) sys->parent;
    DDC* ddc = simulation->ddc;
    if (ddc->lastUpdate == sys->loop)
    { // re-assign map upon domain changes
        assignRestraintMap(sys, parms);
    }

    double xBox = sys->box->h0.xx;
    double yBox = sys->box->h0.yy;
    double zBox = sys->box->h0.zz;
    double xBoxHalf = xBox / 2.0;
    double yBoxHalf = yBox / 2.0;
    double zBoxHalf = zBox / 2.0;

    //    scaleRestraintByBoxChange(parms);

    PUSH_RANGE("restraint", 4);
    int nRestraint = parms->nRestraint;

    double eRestraint = 0.0;
    for (int r = 0; r < nRestraint; r++)
    {
        RESTRAINTPARMS *restraintparms = parms->restraintList[r];

        int ii = parms->restraintMap[r];

        if (ii >= 0)
        {
            //           printf("%d specie=%s gid=%llu ii=%d x=%f y=%f z=%f r=%d x0=%f y0=%f z0=%f\n", sys->loop, state->species[ii]->name, 
            //                state->label[ii], ii, state->rx[ii], state->ry[ii], state->rz[ii],
            //                r, restraintparms->x0, restraintparms->y0, restraintparms->z0);

            //Get the X0, Y0, Z0


            double x0 = restraintparms->x0 * xBox;
            double y0 = restraintparms->y0 * yBox;
            double z0 = restraintparms->z0 * zBox;
            // Assume we are using the Box with the origin at (0, 0, 0)
            if (parms->origin == 0)
            {
                x0 -= xBoxHalf;
                y0 -= yBoxHalf;
                z0 -= zBoxHalf;
            }

            double kb = restraintparms->kb;

            double xDelta = state->rx[ii] - x0;
            double yDelta = state->ry[ii] - y0;
            double zDelta = state->rz[ii] - z0;

            if (restraintparms->fcx > 0 && fabs(xDelta) > xBoxHalf)
            {
                nearestImage(&xDelta, &yDelta, &zDelta);
            }

            if (restraintparms->fcy > 0 && fabs(yDelta) > yBoxHalf)
            {
                nearestImage(&xDelta, &yDelta, &zDelta);
            }

            if (restraintparms->fcz > 0 && fabs(zDelta) > zBoxHalf)
            {
                nearestImage(&xDelta, &yDelta, &zDelta);
            }

            double cxDelta = restraintparms->fcx*xDelta;
            double cyDelta = restraintparms->fcy*yDelta;
            double czDelta = restraintparms->fcz*zDelta;

            eRestraint = eRestraint + kb * (cxDelta * xDelta + cyDelta * yDelta + czDelta * zDelta);

            double kforce = -2 * kb;
            double fxDelta = kforce*cxDelta;
            double fyDelta = kforce*cyDelta;
            double fzDelta = kforce*czDelta;
            //printf("r %i ii %i kb %f e %f fx %f %f %f xd %f %f %f  fcx %i %i %i\n", r,ii,kb, kb*(cxDelta*xDelta+cyDelta*yDelta+czDelta*zDelta), fxDelta, fyDelta, fzDelta, xDelta, yDelta, zDelta, restraintparms->fcx, restraintparms->fcy, restraintparms->fcz);
            //printf("r %i ii %i e %.12f f %.12f %.12f %.12f c %.12f %.12f %.12f kforce %f !\n\n\n", r, ii,kb*(cxDelta*xDelta+cyDelta*yDelta+czDelta*zDelta), fxDelta, fyDelta, fzDelta,cxDelta, cyDelta, czDelta, kforce);
            // printf("rx %.12f %.12f %.12f restr %.12f %.12f %.12f x %.12f %.12f %.12f \n", state->rx[ii], state->ry[ii], state->rz[ii],restraintparms->x0, restraintparms->y0,restraintparms->z0, xDelta,yDelta,zDelta);

            //exit(0);
            state->fx[ii] += fxDelta;
            state->fy[ii] += fyDelta;
            state->fz[ii] += fzDelta;

            // Pressure
            e->virial.xx += fxDelta*cxDelta;
            e->virial.xy += fxDelta*cyDelta;
            e->virial.xz += fxDelta*czDelta;
            e->virial.yy += fyDelta*cyDelta;
            e->virial.yz += fyDelta*czDelta;
            e->virial.zz += fzDelta*czDelta;
        }
    }
    POP_RANGE();
    e->eion += eRestraint;


}
