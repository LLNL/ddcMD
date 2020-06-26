#include "doubleMirror.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include "state.h"
#include "system.h" 
#include "three_algebra.h" 
#include "object.h" 
#include "preduce.h"
#include "ddcMalloc.h"
#include "units.h"
extern FILE* ddcfile;
#include "external.h"
#include "mpiUtils.h"
#include "format.h"
typedef struct DoubleMirror_parms_st
{
   THREE_VECTOR point1;
   THREE_VECTOR point2;
   THREE_VECTOR normal1;
   THREE_VECTOR normal2;
   THREE_VECTOR dP1;
   THREE_VECTOR dP2;
   THREE_VECTOR localdP1;
   THREE_VECTOR localdP2;
   double dPZeroTime; 
   double       v1; //velocity of mirror 1 (in normal direction)
   double       v2; //velocity of mirror 2 	
   double	time_initial; //initial time
   int outputRate;
   FILE *file; 
} DOUBLE_MIRROR_PARMS;

void DoubleMirror_write_dynamics(GROUP *g, FILE *file) 
{
   DOUBLE_MIRROR_PARMS *p=g->parm ;
   double cl=units_convert(1.0,NULL,"l");
   double cv=units_convert(1.0,NULL,"l/t");
   fprintf(file, "%s %s {\n",g->name,g->objclass);
   fprintf(file, "  type = DOUBLE_MIRROR;\n");
   fprintf(file, "  point1 = %f %f %f;\n",  cl*p->point1.x, cl*p->point1.y, cl*p->point1.z);
   fprintf(file, "  point2 = %f %f %f;\n",  cl*p->point2.x, cl*p->point2.y, cl*p->point2.z);
   fprintf(file, "  normal1 = %f %f %f;\n", p->normal1.x, p->normal1.y, p->normal1.z),
   fprintf(file, "  normal2 = %f %f %f;\n", p->normal2.x, p->normal2.y, p->normal2.z),
   fprintf(file, "  v1 = %f;\n",cv*p->v1);
   fprintf(file, "  v2 = %f;\n",cv*p->v2);
   fprintf(file, "}\n");
} 

void doubleMirror_Update(GROUP *g, int mode, void* stateVoid, double time_in, double dt_in)
{
   DOUBLE_MIRROR_PARMS* parms = (DOUBLE_MIRROR_PARMS*) g->parm; //load mirror parameters 
   	parms->point1.x += parms->v1*parms->normal1.x*dt_in; 
   	parms->point1.y += parms->v1*parms->normal1.y*dt_in; 
   	parms->point1.z += parms->v1*parms->normal1.z*dt_in; 

   	parms->point2.x += parms->v2*parms->normal2.x*dt_in; 
   	parms->point2.y += parms->v2*parms->normal2.y*dt_in; 
   	parms->point2.z += parms->v2*parms->normal2.z*dt_in; 
      
   	backInBox_fast(&parms->point1.x, &parms->point1.y, &parms->point1.z); //put mirrors back in box
   	backInBox_fast(&parms->point2.x, &parms->point2.y, &parms->point2.z);
      //DO I NEED TO DO ANYTHING ABOUT FRONT/BACK TIMESTEPS? 
}
void doubleMirror_Update2(GROUP *g, int mode, void* stateVoid, double time_in, double dt_in)
{
   DOUBLE_MIRROR_PARMS *p=g->parm ;
   SIGNED64  loop = system_getLoop(NULL) ;
   if (TEST0(loop,p->outputRate) )
   {
      double dT = time_in-p->dPZeroTime; 
      THREE_VECTOR local[2],global[2];
      local[0] = p->localdP1; 
      local[1] = p->localdP2; 
      MPI_Reduce(local, global, 6, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
      if (getRank(0) == 0 && dT > 0.0) 
      {
         p->dP1 = global[0];
         p->dP2 = global[1];
         double cF=units_convert(1.0,NULL,"amu*Ang*fs^(-2)");
         THREE_VECTOR dP1 = p->dP1; 
         THREE_VECTOR dP2 = p->dP2; 
         VSCALE(dP1,cF/dT); 
         VSCALE(dP2,cF/dT); 

         fprintf(p->file,loopFormat(),loop); 
         fprintf(p->file, " %10.5f",time_in); 
         fprintf(p->file, " %15.8e %15.8e %15.8e",dP1.x,dP1.y,dP1.z); 
         fprintf(p->file, " %15.8e %15.8e %15.8e\n",dP2.x,dP2.y,dP2.z); 
         fflush(p->file); 
      }
      p->localdP1 = vzero; 
      p->localdP2 = vzero; 
      p->dPZeroTime=time_in; 
   }
}
void doubleMirror_velocityUpdate(int mode, int k, GROUP *g, STATE *state, double time, double dt)
{
   double *rx = state->rx; 
   double *ry = state->ry; 
   double *rz = state->rz; 
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 

   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;
   double a = dt/mass; 


   vx[k] += a*fx[k] ; //calculate new velocity of particle (before impact)
   vy[k] += a*fy[k] ;
   vz[k] += a*fz[k] ;


   DOUBLE_MIRROR_PARMS* parms = (DOUBLE_MIRROR_PARMS*) g->parm; //load mirror parameters 

   THREE_VECTOR point1_new = parms->point1; //current position of mirror 1
   THREE_VECTOR point2_new = parms->point2; //current position of mirror 2

   THREE_VECTOR r1; //distance to mirror 1

   r1.x = rx[k] - point1_new.x; 
   r1.y = ry[k] - point1_new.y;
   r1.z = rz[k] - point1_new.z; 
   nearestImage_fast(&r1.x, &r1.y, &r1.z);

   THREE_VECTOR r2; //distance to mirror 2
   r2.x = rx[k] - point2_new.x; 
   r2.y = ry[k] - point2_new.y; 
   r2.z = rz[k] - point2_new.z; 
   nearestImage_fast(&r2.x, &r2.y, &r2.z);

   int mirrorFlag = 1; 
   THREE_VECTOR normal = parms->normal1;
   THREE_VECTOR rr = r1;
   double       v_mirror = parms->v1;
   if ( fabs(DOT(r1, parms->normal1)) > fabs(DOT(r2, parms->normal2)) )
   {
      //point = parms->point2;
      mirrorFlag = 2; 
      normal = parms->normal2;
      rr = r2;
      v_mirror = parms->v2;
   }

   double dot = DOT(rr, normal);
   if (dot > 0) return; // particle on good side of plane

   double vParallel = vx[k]*normal.x + vy[k]*normal.y + vz[k]*normal.z;
   if ((vParallel-v_mirror) > 0) return;       // particle travelling faster than the mirror
   if (mirrorFlag==1) VSVOP(parms->localdP1,+=,2*mass*(v_mirror-vParallel),*,normal); 
   if (mirrorFlag==2) VSVOP(parms->localdP2,+=,2*mass*(v_mirror-vParallel),*,normal); 
   vx[k] += 2*(v_mirror-vParallel)*normal.x;
   vy[k] += 2*(v_mirror-vParallel)*normal.y;
   vz[k] += 2*(v_mirror-vParallel)*normal.z;
}
/*
   {
   MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   */



/*
   void doubleMirror_velocityUpdate1(int mode, int k, GROUP *g,SPECIES *species,THREE_VECTOR V, THREE_VECTOR Ft1, THREE_VECTOR Rt1, double time, double dt)
   {
   double mass = ((ATOMTYPE_PARMS *) (species->parm))->mass;
//VSVOP(Vt1,+=,dt/mass,*,Ft1);
//printf("F= %f %f %f",F.x,F.y,F.z);

double a = dt/mass; 

THREE_VECTOR V2t=V;
V.x += a*Ft1.x ; //calculate new velocity of particle (before impact)
V.y += a*Ft1.y ;
V.z += a*Ft1.z ;
DOUBLE_MIRROR_PARMS* parms = (DOUBLE_MIRROR_PARMS*) g->parm; //load mirror parameters 

THREE_VECTOR point1_new = parms->point1; //current position of mirror 1
THREE_VECTOR point2_new = parms->point2; //current position of mirror 2

point1_new.x += (time-parms->time_initial)*parms->v1*parms->normal1.x; 
point1_new.y += (time-parms->time_initial)*parms->v1*parms->normal1.y; 
point1_new.z += (time-parms->time_initial)*parms->v1*parms->normal1.z;

point2_new.x += (time-parms->time_initial)*parms->v2*parms->normal2.x; 
point2_new.y += (time-parms->time_initial)*parms->v2*parms->normal2.y; 
point2_new.z += (time-parms->time_initial)*parms->v2*parms->normal2.z; 


backInBox_fast(&point1_new.x, &point1_new.y, &point1_new.z); //put mirrors back in box
backInBox_fast(&point2_new.x, &point2_new.y, &point2_new.z);

//parms->point1c = point1_new;
//parms->point2c = point2_new;
parms->time_c  = time;
THREE_VECTOR r1; //distance to mirror 1

r1.x = Rt1.x - point1_new.x; 
r1.y = Rt1.y - point1_new.y;
r1.z = Rt1.z - point1_new.z; 
nearestImage_fast(&r1.x, &r1.y, &r1.z);

THREE_VECTOR r2; //distance to mirror 2
r2.x = Rt1.x - point2_new.x; 
r2.y = Rt1.y - point2_new.y;
r2.z = Rt1.z - point2_new.z; 
nearestImage_fast(&r2.x, &r2.y, &r2.z);

THREE_VECTOR normal = parms->normal1;
THREE_VECTOR rr = r1;
double       v_mirror = parms->v1;
if ( fabs(DOT(r1, parms->normal1)) > fabs(DOT(r2, parms->normal2)) )
{
//point = parms->point2;
normal = parms->normal2;
rr = r2;
v_mirror = parms->v2;
}

double dot = DOT(rr, normal);
if (dot > 0) return; // particle on good side of plane

double vParallel = V2t.x*normal.x + V2t.y*normal.y + V2t.z*normal.z;
if ((vParallel-v_mirror) > 0) return;       // particle travelling faster than the mirror
V2t.x += 2*(v_mirror-vParallel)*normal.x;
V2t.y += 2*(v_mirror-vParallel)*normal.y;
V2t.z += 2*(v_mirror-vParallel)*normal.z;
return V;
}
*/

void doubleMirror_parms(GROUP *gp)
{
   DOUBLE_MIRROR_PARMS* parms = ddcMalloc(sizeof(DOUBLE_MIRROR_PARMS));

   OBJECT* obj = (OBJECT*) gp;

   object_get(obj, "point1",  &parms->point1,  WITH_UNITS, 3, "0 0 -1", "l", NULL);
   object_get(obj, "point2",  &parms->point2,  WITH_UNITS, 3, "0 0 1",  "l", NULL);
   object_get(obj, "normal1", &parms->normal1, DOUBLE,     3, "0 0 1",     "l", NULL);
   object_get(obj, "normal2", &parms->normal2, DOUBLE,     3, "0 0 -1",    "l", NULL);
   object_get(obj, "v1", &parms->v1, WITH_UNITS,     1, "0.0",     "l/t", NULL);
   object_get(obj, "v2", &parms->v2, WITH_UNITS,     1, "0.0",     "l/t", NULL);
   object_get(obj, "outputRate", &parms->outputRate, INT,  1, "0");

   parms->time_initial=system_getTime((SYSTEM *) gp->parent);
   parms->dP1 = vzero; 
   parms->dP2 = vzero; 
   parms->localdP1 = vzero; 
   parms->localdP2 = vzero; 
   parms->dPZeroTime = parms->time_initial; 
   char filename[256]; 
   assert(strlen(gp->name) < 240); 
   sprintf(filename,"%s_dPdt.data",gp->name); 
   if (getRank(0) ==0 ) 
   {  
      parms->file = fopen(filename,"a");
      fprintf(parms->file,"%-10s %-10s %-45s %-45s\n","#loop"," time (fs)", "            dP1/dt (amu*Ang/fs^2)","            dP2/dt (amu*Ang/fs^2)"); 
      fflush(parms->file); 
      printf("filename=%s\n",filename);fflush(stdout); 
   }

   double tmp = 1.0/sqrt(VSQ(parms->normal1));
   VSCALE(parms->normal1, tmp);
   tmp = 1.0/sqrt(VSQ(parms->normal2));
   VSCALE(parms->normal2, tmp);

   gp->itype = DOUBLE_MIRROR;
   gp->parm = parms;
   gp->write_dynamics = DoubleMirror_write_dynamics;
   gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))doubleMirror_velocityUpdate; 
   gp->Update= (void (*)(GROUP *, int, void *, double, double))doubleMirror_Update; 
   gp->Update2= (void (*)(GROUP *, int, void *, double, double))doubleMirror_Update2; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
