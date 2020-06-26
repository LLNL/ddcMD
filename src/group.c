#include "group.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "ddcMalloc.h"
#include "berendsen.h"
#include "doubleMirror.h"
#include "mpiUtils.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
static GROUP **grouplist = NULL;
static int ngroups = 0;
static char *emptyString = "";
static void badGroupTypeError(const char* name, const char* type);

void free_parms(GROUP *);
void langevin_parms(GROUP *);
void relativisticLangevin_parms(GROUP *);
void TNBurn_parms(GROUP *);
void extforce_parms(GROUP *);
void shear_parms(GROUP *gp);
void shwall_parms(GROUP *gp);
void fixedVelocity_parms(GROUP *);
void relativistic_parms(GROUP *);
void frozen_parms(GROUP *);
void ionization_parms(GROUP *gp){};  //Hack to avoid compiling the non functional ionization code
void radiation_parms(GROUP *);
void piston_parms(GROUP *);
void quench_parms(GROUP*);
void unionGroup_parms(GROUP*);
char *group_parse (struct group_st *g, char *string, int i){ return string;}
char *group_write (struct group_st *g, int i){return emptyString;}
unsigned  group_bread (unsigned char* buf, struct group_st* g, int i){return 0;}
unsigned  group_bwrite (unsigned char* buf, struct group_st *g, int i){return 0;}
void group_Update(GROUP *g, int mode, void *state, double time_in, double dt_in) {  }
// (void (*) (GROUP*, int, void *, double, double))
GROUP *group_init(void *parent,char *name)
{
	GROUP*  gp = (GROUP *) object_initialize(name, "GROUP", sizeof(GROUP));
	gp->parent = parent; 

	char *type;
	object_get((OBJECT *) gp, "type", &type, STRING, 1, "");
	object_get((OBJECT *) gp, "volume_per_atom", &gp->energyInfo.vol_per_atom, DOUBLE, 1, "-1");
	gp->type = strdup(type);
	gp->index = ngroups; 
	gp->defaultValue = NULL;
	gp->checkValue = NULL;
	gp->parse = NULL;
	gp->write = NULL;
	gp->bwrite = NULL;
	gp->write_dynamics = NULL;
	gp->bread = NULL; 
	gp->Update =group_Update;
   gp->Update1=group_Update; 
   gp->Update2=group_Update; 
	gp->start = group_Update;
	gp->velocityUpdate = NULL;
	gp->velocityUpdateKernel = NULL;
	gp->nMember = 0;
	gp->useDefault = 0;
	gp->parm = NULL;
	gp->energyBath = 0.0;
	gp->energyBathFunction = NULL;
   gp->energyInfo.virialCorrection = NULL; 
	grouplist = ddcRealloc(grouplist, sizeof(GROUP *)*(ngroups + 2));
	grouplist[ngroups++] = gp;
	grouplist[ngroups] = NULL;
	if (strcmp(type, "BERENDSEN") == 0) berendsen_parms(gp); 
	if (strcmp(type, "EXTFORCE") == 0) extforce_parms(gp);
	if (strcmp(type, "FIXEDVELOCITY") == 0) fixedVelocity_parms(gp); 
	if (strcmp(type, "RELATIVISTIC") == 0) relativistic_parms(gp); 
	if (strcmp(type, "FREE") == 0) free_parms(gp);
	if (strcmp(type, "FROZEN") == 0) frozen_parms(gp);
	if (strcmp(type, "IONIZATION") == 0) ionization_parms(gp); 
	if (strcmp(type, "LANGEVIN") == 0) langevin_parms(gp);
	if (strcmp(type, "RELATIVISTICLANGEVIN") == 0) relativisticLangevin_parms(gp);
	if (strcmp(type, "TNBURN") == 0) TNBurn_parms(gp);
	if (strcmp(type, "PISTON") == 0) piston_parms(gp);
	if (strcmp(type, "QUENCH") == 0) quench_parms(gp); 
	if (strcmp(type, "RADIATION") == 0) radiation_parms(gp); 
	if (strcmp(type, "SHEAR") == 0) shear_parms(gp);
	if (strcmp(type, "SHWALL") == 0) shwall_parms(gp); 
	if (strcmp(type, "DOUBLE_MIRROR") == 0) doubleMirror_parms(gp); 
	if (strcmp(type, "UNIONGROUP") == 0) unionGroup_parms(gp); 
//  Next block is old dinosar tracks.  We keep it around in case we ever
//  want these thermostats.
/* 	if (strcmp(type, "ANDERSON") == 0) */
/* 	{ */
/* 		ANDERSEN_PARMS *parms; */
/* 		parms = ddcCalloc(1, sizeof(ANDERSEN_PARMS)); */
/* 		gp->itype = ANDERSON; */
/* 		nelements = object_get((OBJECT *) gp, "temperature", &temp_str, STRING, 1, "300.0"); */
/* 		nelements = object_get((OBJECT *) gp, "dTdt", &parms->dTdt, DOUBLE, 1, "0.0"); */
/* 		gp->parm = (void *)parms; */
/* 	} */
/* 	if (strcmp(type, "NOSE") == 0) */
/* 	{ */
/* 		NOSETYPE_PARMS *parms; */
/* 		parms = ddcCalloc(1, sizeof(NOSETYPE_PARMS)); */
/* 		gp->itype = NOSE; */
/* 		nelements = object_get((OBJECT *) gp, "method", &parms->method, STRING, 1, "KINETIC"); */
/* 		nelements = object_get((OBJECT *) gp, "tau", &parms->tau, DOUBLE, 1, "1.0"); */
/* 		nelements = object_get((OBJECT *) gp, "zeta", &parms->zeta, DOUBLE, 1, "0.0"); */
/* 		nelements = object_get((OBJECT *) gp, "temperature", &temp_str, STRING, 1, "300.0"); */
/* 		nelements = object_get((OBJECT *) gp, "dTdt", &parms->dTdt, DOUBLE, 1, "0.0"); */
/* 		nelements = object_get((OBJECT *) gp, "v", &parms->v, DOUBLE, 3, "0.0"); */
/* 		nelements = object_get((OBJECT *) gp, "zeta1", &zeta1, DOUBLE, 1, "0.0"); */
/* 		if (strcmp(parms->method, "KINETIC") == 0) parms->imethod = KINETIC; */
/* 		if (strcmp(parms->method, "VELOCITY") == 0) parms->imethod = VELOCITY; */
/* 		sscanf(temp_str, "%lf%s", &parms->Target_temperature, units); */
/* 		if (strcmp(units, "K") == 0) */
/* 		{ */
/* 			parms->Target_temperature /= (11605.0); */
/* 			parms->dTdt /= (11605.0); */
/* 		} */
/* 		gp->parm = (void *)parms; */
/* 	} */
	if (gp->Update == NULL || gp->velocityUpdate == NULL) badGroupTypeError(name, type);

	return gp;
}

// more dinosaur tracks
/* void nose_write_dynamics(GROUP* gp, FILE* file) */
/* { */
/*    fprintf(file, "%s %s { zeta = %e ;}\n", */
/* 	   gp->name, gp->objclass, ((NOSETYPE_PARMS *) gp->parm)->zeta); */
/* } */

void group_get(GROUP*group, int get, void **ptr)
{
	int i, max;
	switch (get)
	{
	case GROUPTYPE:
		*((char **)ptr) = group->type;
		return;
	case GROUPITYPE:
		*((int *)ptr) = group->itype;
		return;
	case NGROUPS:
		*(int *)ptr = ngroups;
		return;
	case GROUPLIST:
		*((GROUP ***) ptr) = grouplist;
		return;
	case GROUPMAXNAMELENGTH:
		i = max = 0;
		while (grouplist[i] != NULL)
		{
			max = MAX(max, (int)strlen(grouplist[i]->name));
			i++;
		}
		*(int *)ptr = max;
		return;
	case GROUPMAXWRITELENGTH:
		i = max = 0;
		while (grouplist[i] != NULL)
		{
			if (grouplist[i]->write != NULL) max = MAX(max, (int)strlen(grouplist[i]->write(grouplist[i], -1)));
			i++;
		}
		*(int *)ptr = max;
		return;
	case GROUPMAXBINARYWRITELENGTH:
		i = max = 0;
		while (grouplist[i] != NULL)
		{
			if (grouplist[i]->bwrite != NULL) max = MAX(max, (int)grouplist[i]->bwrite(NULL,grouplist[i], -1));
			i++;
		}
		*(int *)ptr = max;
		return;
	default:
		break;
	}
}

int group_put(GROUP*group, int put, void **ptr)
{
	switch (put)
	{
	default:
		return 0;
	}
}

GROUP *group_find(GROUP ** group, char *name)
{
	int i;
	i = 0;
	if (group == NULL) group = grouplist;
	while (group[i] != NULL)
	{
		if (strcmp(group[i]->name, name) == 0) return group[i];
		i++;
	}
	if (i==1) return group[0];
	return NULL; 
}

GROUP *group_by_index(GROUP ** group, int index)
{
	int i;
	i = 0;
	if (group == NULL) group = grouplist;
	while (group[i] != NULL)
	{
		if (group[i]->index == index) return group[i];
		i++;
	}
	return NULL;
}

void badGroupTypeError(const char* name, const char* type)
{
   if (getRank(0) != 0) return;

   if (strlen(type) == 0)
      printf("group_init:  It appears you forgot to specify the type for the %s group\n", name);
   else
      printf("group_init: Unrecognized group type %s for the %s group\n", type, name);
   abortAll(6);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
