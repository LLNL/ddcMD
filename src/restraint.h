#ifndef RESTRAINT_H
#define RESTRAINT_H

#include "gid.h"
#include "system.h"
#include "cudaUtils.h"

typedef struct restraintparms_st
{
    char *name;
    char *objclass;
    char *value;
    char *type; /* model */
    void *parent;

    gid_type gid;

    int atomI;
    int func;

    int fcx;
    int fcy;
    int fcz;

    double x0;
    double y0;
    double z0;

    double kb;

} RESTRAINTPARMS;

typedef struct restraint_gpu_parms_st
{
    unsigned *restraintMap;

    gid_type gid;

    int atomI;
    int func;

    int *fcx;
    int *fcy;
    int *fcz;

    double *x0;
    double *y0;
    double *z0;

    double *kb;

} RESTRAINTGPUPARMS;

typedef struct restraint_st
{
    char *name;
    char *objclass;
    char *value;
    char *type; /* model */
    void *parent;

    int origin;
    int nRestraint;
    RESTRAINTPARMS** restraintList;
    gid_type *restGidList; // A gid list helps binary search

    unsigned *restraintMap;

    RESTRAINTGPUPARMS *gpu_parms;

} RESTRAINTLIST;

#ifdef __cplusplus
extern "C"
{
#endif


GPUFUNC(void restraintGPU_parms(RESTRAINTLIST * restraint))
GPUFUNC(void restraintGPU(SYSTEM *sys, RESTRAINTLIST *parms, ETYPE *e))


void assignRestraintMap(SYSTEM *sys, RESTRAINTLIST *restraint);
#ifdef __cplusplus
}
#endif


#endif /* RESTRAINT_H */

