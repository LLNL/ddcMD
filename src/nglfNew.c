#include "nglfNew.h"
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
#include "constraint.h"
void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
void solveConstraintGoup(CONSTRAINTNEW *constraint, double dt, STATE *state, THREE_SMATRIX *sion, int location);
static THREE_SMATRIX velocityConstraint(double dt, STATE *state, NGLFNEW_PARMS *parms,  int location);
static double frontFunc(double dt, double d2, THREE_VECTOR rab,THREE_VECTOR vab);
static double backFunc(double dt, double d2, THREE_VECTOR rab,THREE_VECTOR vab);

//static AUXNEIGHBOR *auxNeighbor;

static void adjustPosn(STATE *state, BOX_STRUCT* box) 
{
   THREE_MATRIX hfac; 
   box_get(box,HFAC,(void *)&hfac);
   for (int kk = 0; kk < state->nlocal; kk++)
   {
      THREE_VECTOR old={state->rx[kk],state->ry[kk],state->rz[kk]}; 
      THREE_VECTOR new = matrix_vector(hfac,old); 
      state->rx[kk]=new.x; state->ry[kk]=new.y; state->rz[kk]=new.z;
   }
}
static void  changeVolume(STATE *state, BOX_STRUCT *box, THREE_SMATRIX *sion, double P0,double beta, double tau , double dt)
{
   THREE_MATRIX  lambda=mzero; 
   double btt = beta*dt/tau; 
   double Pxx,Pyy,Pzz; 
   Pxx = Pyy = -0.5*(sion->xx+sion->yy); 
   Pzz =        -sion->zz; 

   lambda.xx = cbrt(1.0+(Pxx-P0)*btt); 
   lambda.yy = cbrt(1.0+(Pyy-P0)*btt); 
   lambda.zz = cbrt(1.0+(Pzz-P0)*btt); 

   THREE_MATRIX h0=box_get_h(box); 
   THREE_MATRIX h = matrix_matrix(lambda,h0) ; 
   box_put(box,HO,&h); 
   adjustPosn(state,box); 
      
}
NGLFNEW_PARMS *nglfNew_parms(INTEGRATOR*integrator) 
{
    NGLFNEW_PARMS *parms = ddcMalloc(sizeof (NGLFNEW_PARMS));
    object_get((OBJECT*)integrator,"P0",&parms->P0,WITH_UNITS,1,"0.0","pressure",NULL); 
    object_get((OBJECT*)integrator,"beta",&parms->beta,WITH_UNITS,1,"0.0","1/pressure",NULL); 
    object_get((OBJECT*)integrator,"tauBarostat",&parms->tauBarostat,WITH_UNITS,1,"0.0","t",NULL); 
//  double rmax = 3.0;
//  auxNeighbor = auxNeighbor_request(rmax);
//  SIMULATE *simulate=(SIMULATE*)integrator->parent;
//  SYSTEM *system=simulate->system;
    return parms; 
}
static double frontFunc(double dt, double d2, THREE_VECTOR rab,THREE_VECTOR vab) 
{
            THREE_VECTOR pab = rab;
            VSVOP(pab, +=, dt,*,vab); 
            double pab2 = VSQ(pab); 
            double rvab =  (pab2-d2)/(2*dt); 
//          p2-d2 = (r + v*dt)^2 - d2  = 0.5*(r2-d2)/dt + r.v +0.5*dt*v2 
            return rvab; 
}
static double backFunc(double dt, double d2, THREE_VECTOR rab,THREE_VECTOR vab) 
{
            double rvab = DOT(rab,vab); 
            return rvab; 
}

void solveConstraintGroup(CONSTRAINTNEW *constraint, double dt, STATE* state, THREE_SMATRIX *sion, int location)
{

   double ( *func) (double dt, double d2, THREE_VECTOR rab,THREE_VECTOR vab);
   if (location == FRONT_TIMESTEP) func = frontFunc; 
   if (location == BACK_TIMESTEP ) func = backFunc; 
   const double tol=1.0e-8; //unit less
   const int maxit=500;

   int numAtm = constraint->numAtom; // total number of atom
   double   rMass[numAtm]; 
   THREE_VECTOR r[numAtm],v[numAtm];
   for(int j = 0; j < numAtm; ++j)
   {
      int index=constraint->atomList[j];
      rMass[j]=1.0/((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
      VSET(r[j],state->rx[index],state->ry[index],state->rz[index]);
      VSET(v[j],state->vx[index],state->vy[index],state->vz[index]);
   }
   int numPair = constraint->numPair; 
   THREE_VECTOR rab[numPair];
   double gamma[numPair]; 
   for (int ab = 0; ab < numPair; ++ab) 
   {
      CONSPAIR* consPair=constraint->consPair[ab];
      int a=consPair->atomI;
      int b=consPair->atomJ;
      VOP2(rab[ab],=,r[a],-,r[b]); 
      nearestImage(&rab[ab].x, &rab[ab].y, &rab[ab].z);
   }
   double errMax = 0.0;    
   int it=0;
   for (int ab = 0; ab < numPair; ++ab) gamma[ab]=0.0; 
   for (;it<maxit;it++) //start the iterative loop
   {
      errMax = 0.0;    
      for (int ab = 0; ab < numPair; ++ab) 
      {
         CONSPAIR* consPair=constraint->consPair[ab];
         int a=consPair->atomI;
         int b=consPair->atomJ;
         double d2=consPair->d2;
         THREE_VECTOR vab; 
         VOP2(vab,=,v[a],-,v[b]); 
         double rma=rMass[a];
         double rmb=rMass[b];                       

         double rvab = func(dt,d2, rab[ab], vab)/d2 ; 

         double gab=-rvab/(rma+rmb);   //units: mass/time)
         double err = fabs(rvab*dt);
         if (err > errMax ) errMax = err;  
         VSVOP(v[a],+=,(rma*gab),*,rab[ab]); 
         VSVOP(v[b],-=,(rmb*gab),*,rab[ab]); 
         gamma[ab] += gab; 
         //if (err < tol) printf("%d: group=%d ab= %d numPair=%d it=%d %e %e\n",getRank(0),i,ab,numPair,it,err, gab); 
      } 
      if (errMax < tol ) 
      {
         break;   
      }
   } //end of iterative loop
   if(it == maxit) printf("%d SolveConstraintGroup: too many contraint iterations.\n",getRank(0));
   *sion = szero; 
   for (int ab = 0; ab < numPair; ++ab) 
   {
      THREE_VECTOR fab;
      VSVOP(fab,=,(-2.0/dt *gamma[ab]),*,rab[ab]);  //fab = gamma/(dt/2)*rab. delta moment v*m = f *t ;  In this case t = dt/2
      sion->xx+=fab.x*rab[ab].x; 
      sion->yy+=fab.y*rab[ab].y; 
      sion->zz+=fab.z*rab[ab].z; 
      sion->xy+=fab.x*rab[ab].y; 
      sion->xz+=fab.x*rab[ab].z; 
      sion->yz+=fab.y*rab[ab].z; 
   }
   for(int j = 0; j < numAtm; ++j) // Store the new values
   {
      int index=constraint->atomList[j];
      state->vx[index]=v[j].x;
      state->vy[index]=v[j].y;
      state->vz[index]=v[j].z;   
   }            
}

static THREE_SMATRIX velocityConstraint(double dt, STATE *state, NGLFNEW_PARMS *parms, int location)
{
   THREE_SMATRIX sion=szero; 
   for (int i=0; i < parms->nConstraintGroup; i++) solveConstraintGroup(parms->constraints[i], dt, state, &sion, location);
   return sion; 
}
void positionUpate(double dt, STATE *state)
{
   //REF *v0 = vaf_v0(); 
   for (int kk = 0; kk < state->nlocal; kk++)
   {
      THREE_VECTOR delta; 
      state->rx[kk] += delta.x = dt*state->vx[kk];
      state->ry[kk] += delta.y = dt*state->vy[kk];
      state->rz[kk] += delta.z = dt*state->vz[kk];
//      if (v0 != NULL )  VOP1(v0[kk].r,+=,delta);   // for GPU code you could assume v0 == NULL;  can there VOP1 and delta are not needed. Shiv talk to me about this.  Jim 
   }
}
void nglfNew(DDC*ddc, SIMULATE*simulate, NGLFNEW_PARMS* parms ) 
{
   double dt = simulate->dt;
   double time = simulate->time;
   SYSTEM* sys = simulate->system;
   STATE* state = sys->collection->state;
   //scalePositionsByBoxChange_sk(sys->box, time, state->rx, state->ry, state->rz, state->nlocal);
   changeVolume(state, sys->box, &(sys->energyInfo.sion), parms->P0, parms->beta, parms->tauBarostat , dt);
   for (int ii=0;ii<sys->ngroup;ii++) sys->group[ii]->velocityUpdateKernel(FRONT_TIMESTEP,sys->group[ii],state,time,0.5*time); 
   THREE_SMATRIX sionFront = velocityConstraint(dt, state, parms, FRONT_TIMESTEP);

   positionUpate(dt, state);

   //scalePositionsByBoxChange_sk(sys->box, time+dt, state->rx, state->ry, state->rz, state->nlocal);
   for (int kk = 0; kk < state->nlocal; kk++) backInBox_fast(state->rx + kk, state->ry + kk, state->rz + kk);

   ddc->update = 0;
   time += dt; // positions, box (volume, h0, hinv), and forces at  t = n*dt + dt 
   simulate->time = sys->time = time;
   simulate->loop++;
   sys->loop = simulate->loop;

   for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update1(sys->group[kk], -1, state, time, 0.5 * dt);

   if (ddcenergy(ddc, sys, 0) != 0) return;

   for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update(sys->group[kk], BACK_TIMESTEP, state, time, 0.5 * dt);
   for (int ii=0;ii<sys->ngroup;ii++) sys->group[ii]->velocityUpdateKernel(BACK_TIMESTEP,sys->group[ii],state,time,0.5*time); 

   THREE_SMATRIX sionBack = velocityConstraint(dt, state, parms, BACK_TIMESTEP);

   sys->energyInfo.sion.xx+=0.5*(sionFront.xx+sionBack.xx);    //Note sion needs to be  will be accumlated over task and/or threads, but the global value of sion  can be updated until 
   sys->energyInfo.sion.yy+=0.5*(sionFront.yy+sionBack.yy);    //after the volume update step has been completed.  
   sys->energyInfo.sion.zz+=0.5*(sionFront.zz+sionBack.zz);    //Need to think about restart.   With the current implementation there is no constraint correction on the first timestep.  
   sys->energyInfo.sion.xy+=0.5*(sionFront.xy+sionBack.xy);    //so on every first restart there there is a small error as well reproducability is broken. 
   sys->energyInfo.sion.xz+=0.5*(sionFront.xz+sionBack.xz); 
   sys->energyInfo.sion.yz+=0.5*(sionFront.yz+sionBack.yz); 

   kinetic_terms(sys, 1);

   for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update2(sys->group[kk], -1, state, time, 0.5 * dt);
   //eval_energyInfo(sys);
   for (int kk = 0; kk < sys->ngroup; kk++) sys->group[kk]->Update(sys->group[kk], FRONT_TIMESTEP, state, time, 0.5 * dt);

   simulate->time = sys->time = time;
   /*errorCheck(ddc->domain_id, simulate->loop, state, sys->energyInfo, parms, datafile); */
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
