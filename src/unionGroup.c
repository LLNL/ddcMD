//#include "unionGroup.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include "object.h"
#include "group.h"
#include "error.h"
#include "ddcMalloc.h"
/*
 *
 * Groups may also define the following optional functions:
 *    void (*defaultValue)      (struct group_st *g, gid_type label, int i);
 *    int  (*checkValue)         (struct group_st *g, int i);
 *    char*     (*parse)        (struct group_st *g, char *string, int i);
 *    char*    (*write)         (struct group_st *g, int i);
 *    unsigned (*bread)         (unsigned char* buf, struct group_st* g, int i);
 *    unsigned (*bwrite)        (unsigned char* buf, struct group_st *g, int i);
 *    void     (*write_dynamics)(struct group_st *g, FILE *file);
 *    void     (*start)         (struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*Update)        (struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*Update1)(struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*Update2)(struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*velocityUpdate)(int mode, int k, struct group_st *g, STATE *state, double time, double dt);
 *
*/
typedef struct unionGroupParms_str { int ngroups; GROUP **groups;} UNIONGROUP_PARMS; 

void     unionGroup_defaultValue  (GROUP *g, gid_type label, int i)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->defaultValue != NULL) groups[i]->defaultValue(groups[i],label,i); 
}
int      unionGroup_checkValue    (GROUP *g, int i)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   unsigned check=0; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->checkValue != NULL) check = check | groups[i]->checkValue(groups[i],i); 
   return check; 
}
char *unionGroup_parse(GROUP*g, char *field, int i)
{
	UNIONGROUP_PARMS *parms=g->parm ;  
   GROUP **groups = parms->groups; 
   field = strdup(field); 
   char *tok = strtok(field,"|"); 
   int k=0; 
   while(tok != NULL) 
   {
      GROUP *group = groups[k++]; 
      if (group->parse != NULL) group->parse(group,tok,i); 
      tok = strtok(NULL,"|"); 
   }
   free(field); 
	return NULL; 
   
}
char *unionGroup_write(GROUP*g, int i)
{
   static char *line=NULL; 
   static int lineLength=0; 
   UNIONGROUP_PARMS *parms=g->parm ;  
   GROUP **groups = parms->groups; 

   char *string[parms->ngroups]; 
   int cnt=1;  //For terminator
   for (int i=0;i<parms->ngroups;i++) 
   {
      if (groups[i]->write == NULL) string[i] = groups[i]->write(groups[i],i); 
      cnt += strlen(string[i])+1; //extra byte for separater
      assert(index(string[i],'|') != NULL);   // The pipe symbol '|' is not allowed string. 
   }
   if (cnt >  lineLength) 
   {
      line=ddcRealloc(line,cnt); 
      lineLength = cnt; 
   }
   strcpy(line,string[0]); 
   for (int i=1;i<parms->ngroups;i++) 
   {
      strcat(line,"|"); 
      strcat(line,string[i]); 
   }
   return  line; 
}
unsigned unionGroup_bread         (unsigned char* buf, GROUP *g, int i)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   unsigned cnt=0; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->bread != NULL) cnt+=groups[i]->bread(buf+cnt,groups[i],i); 
   return cnt;
}
unsigned unionGroup_bwrite        (unsigned char* buf, GROUP *g, int i)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   unsigned cnt=0; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->bwrite != NULL) cnt+=groups[i]->bwrite(buf+cnt,groups[i],i); 
   return cnt;
}
void     unionGroup_write_dynamics(GROUP *g, FILE *file)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->write_dynamics != NULL) groups[i]->write_dynamics(groups[i],file); 
}
void     unionGroup_Update        (GROUP *g, int mode, void *state, double t, double dt)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->Update != NULL) groups[i]->Update(groups[i],mode,state,t,dt); 
}
void     unionGroup_Update1       (GROUP *g, int mode, void *state, double t, double dt)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->Update1 != NULL) groups[i]->Update1(groups[i],mode,state,t,dt); 
}
void     unionGroup_Update2       (GROUP *g, int mode, void *state, double t, double dt)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->Update2 != NULL) groups[i]->Update2(groups[i],mode,state,t,dt); 
}
void     unionGroup_start         (GROUP *g, int mode, void *state, double t, double dt)
{
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
   for (int i=0;i<parms->ngroups;i++) if (groups[i]->start != NULL) groups[i]->start(groups[i],mode,state,t,dt); 
}
void unionGroup_velocityUpdate(int mode, int k, GROUP *g, STATE *state, double time, double dt)

{	
   UNIONGROUP_PARMS *parms=g->parm;  
   GROUP **groups = parms->groups; 
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
  	vy[k] += a*fy[k] ;
   double v0x = vx[k]; 
   double v0y = vy[k]; 
   double v0z = vz[k]; 
   int aa=1,bb=0; 
   switch(mode)
   {
      case FRONT_TIMESTEP:
      aa = 1;
      bb = 0; 
      break; 
      case BACK_TIMESTEP:
      aa = -1;
      bb = parms->ngroups-1; 
      break; 
   }
   double dvx, dvy, dvz; 
   dvx = dvy = dvz = 0.0; 
   for (int i=0;i<parms->ngroups;i++) 
   {
      int ii = aa*i + bb ; 
      groups[ii]->velocityUpdate(mode,k,groups[ii],state,time,dt); 
      dvx += vx[k] - v0x - a*fx[k]; 
      dvy += vy[k] - v0y - a*fy[k];
      dvz += vy[k] - v0z - a*fz[k];
      vx[k] = v0x; 
      vy[k] = v0y; 
      vz[k] = v0z; 
   }
   vx[k] += dvx + a*fx[k]; 
   vy[k] += dvy + a*fy[k]; 
   vz[k] += dvz + a*fz[k]; 
}

void unionGroup_parms(GROUP *gp)
{
   UNIONGROUP_PARMS  *parms = ddcCalloc(1, sizeof(UNIONGROUP_PARMS));
   gp->itype = UNIONGROUP;
   char **groupNames; 
   parms->ngroups = object_getv((OBJECT *) gp, "groups", (void *)&groupNames, STRING, ABORT_IF_NOT_FOUND);
   parms->groups = ddcMalloc(parms->ngroups*sizeof(GROUP*)); 
   for (int i=0;i<parms->ngroups;i++)
   {
      parms->groups[i] = group_init(gp, groupNames[i]);
   }
   gp->parm = parms;
   gp->defaultValue = unionGroup_defaultValue;
   gp->checkValue = unionGroup_checkValue;
   gp->parse = unionGroup_parse;
   gp->write = unionGroup_write;
   gp->bread = unionGroup_bread;
   gp->bwrite = unionGroup_bwrite;
   gp->write_dynamics = unionGroup_write_dynamics;
   gp->start= unionGroup_start; 
   gp->Update= unionGroup_Update; 
   gp->Update1=unionGroup_Update1; 
   gp->Update2=unionGroup_Update2; 
   gp->velocityUpdate=  (void (*)(int, int, GROUP *,void *,double, double))unionGroup_velocityUpdate; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
