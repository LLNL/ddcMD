#include "ddcenergy.h"

#ifdef BGL
#include "ddcenergy.h"
#include <rts.h>
#endif
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "heap.h"
#include "particle.h"
#include "parityHandler.h"
#include "utilities.h"
#include "collection.h"
#include "neighbor.h"
#include "ptiming.h"
#include "state.h"
#include "auxNeighbor.h"
#include "box.h"
#include "potential.h"
#include "units.h"
#include "loadBalance.h"
#include "simulate.h"
#include "mpiUtils.h"

#define INDEX 0 /* PAIR FUNCTION OPT INDEX  1 */ 
#define MAX(A, B) ((A) > (B) ? (A) : (B))
void makeImageParticles(SYSTEM *sys, double rInfluence); 


LONG64 mgptTempStorageSize(SYSTEM*sys, void*parms);


static int parityCheckingBarrier(MPI_Comm);
static void idleThreadsCheck(SYSTEM*system, MPI_Comm comm);
void fsumX(SYSTEM*sys, THREE_SMATRIX*restrict virial);
void cutoffs(DDC*ddc, SYSTEM*sys)
{
	neighborCutoff(sys->neighbor, -1, (RCUT_TYPE *)NULL);
	ddc->rInfluence = 0.0;
	for (int i=0; i<sys->npotential; ++i)
	{
	   int n;
	   RCUT_TYPE* rcut = sys->potential[i]->getCutoffs(sys, sys->potential[i]->parms, &n);
	   ddc->rInfluence = MAX(ddc->rInfluence, rcut[n].value);
	   neighborCutoff(sys->neighbor, n, rcut);
	}
	ddc->rcut = ddc->rInfluence+sys->neighbor->deltaR; 
}
static LONG64 neighborTempStorageSize(SYSTEM*sys, NBR*nbr)
{
	static RCUT_TYPE *rcut;
	int i;
	LONG64 TempStorageSize; 
	double rlocal2,rremote2,voln; 
	voln = sys->box->volume/(sys->nglobal);
	rlocal2= rremote2 = 0.0; 
	rcut =  nbr->rcut;
	for (i=0;i<nbr->nc;i++) 
	{
		if (rcut[i].mode == RCUT_LOCAL || rcut[i].mode == RCUT_ALL ) 
		{
			rlocal2 = (rcut[i].value)*(rcut[i].value); 
			i++; 
			break; 
		}
	}
	for (;i<nbr->nc;i++) if (rcut[i].mode == RCUT_ALL) {rremote2 = (rcut[i].value)*(rcut[i].value); break; }
	TempStorageSize =       (sys->nlocal)*(2*M_PI/(3.0*voln)*CUBE(sqrt(rlocal2)))*sizeof(RPAIR);
	TempStorageSize += 4.0*(sys->nremote)*(2*M_PI/(3.0*voln)*CUBE(sqrt(rremote2)))*sizeof(RPAIR);
	return  TempStorageSize; 
}
void tempstorage(DDC*ddc, SYSTEM*sys)
{
	int i;
	LONG64 TempStorageSize; 
	double fudge_factor = 1.4; 
	TempStorageSize=0; 
	int needs_heap = 0;
	for (i = 0; i < sys->npotential; i++)
	{
		switch (sys->potential[i]->itype)
		{
		case MGPT:
			//TempStorageSize = MAX(TempStorageSize,mgptTempStorageSize(sys,(MGPT_PARMS *)sys->potential[i]->parms)); 
			TempStorageSize = MAX(TempStorageSize,mgptTempStorageSize(sys,sys->potential[i]->parms)); 
			needs_heap = 1;
			break;
		case EAM:
			needs_heap = 1;
			break;
		case ORDERSH:
			needs_heap = 1;
			break;
		default:
			break;
		}
	}
	
	TempStorageSize += neighborTempStorageSize(sys,sys->neighbor); 
	TempStorageSize *= fudge_factor; 
	TempStorageSize = MAX(32*1024*1024,TempStorageSize); 

	if (needs_heap)
	{
	   if (getRank(0) ==0 && TempStorageSize>32*1024*1024) 
	   {
	      printf("Size of TempStorage in heap: size=%"PRIu64" kb\n",TempStorageSize/1024);
	   }
	   heap_allocate((unsigned int)TempStorageSize); 
	}
}
void zeroParticle(SYSTEM*sys)
{
	double* fx = sys->collection->state->fx;
	double* fy = sys->collection->state->fy;
	double* fz = sys->collection->state->fz;
	double* energy = sys->collection->state->potentialEnergy;
	double* virial = sys->collection->state->virial;
	THREE_SMATRIX *sion  = sys->collection->state->sion;
	if (fx !=NULL) 
	{
		for (unsigned  i = 0; i < sys->nion; i++) 
      {
          fx[i] = fy[i] = fz[i] = 0.0; 
		    energy[i]=virial[i]= 0.0; 
          sion[i] = szero;
      }
	}
}
void zeroEType(ETYPE *e) 
{
	e->eion = e->pion = e->pdiaV = e->e0_at_V = 0.0; 
	e->eion_local = e->pion_local = e->pdiaV = e->e0_at_V = 0.0; 
   e->sion_local=e->sion=e->virial=e->tion=szero;
   e->vcm=e->thermal_flux=vzero; 
	e->rk_local=e->rk=0; 
	e->mass=0.0; 
	e->number=0.0; 
	e->eBath=0.0; 
	e->f = vzero; 
}
void zeroAll(SYSTEM*sys)
{
   sys->nIdleThreads = 0;
   ACCELERATOR *accelerator = accelerator_getAccelerator(NULL);
   // HERE assume if no accelerator is used, the accelerator object will not be created.
   if (accelerator == NULL) zeroParticle(sys); 
   zeroEType(&sys->energyInfo);
   for (int i=0;i<sys->ngroup;i++)   zeroEType(&(sys->group[i]->energyInfo));
   for (int i=0;i<sys->nspecies;i++) zeroEType(&(sys->species[i]->energyInfo));
}

int ddcenergy(DDC*ddc, SYSTEM*sys, int e_eval_flag)
{
	STATE *state=sys->collection->state;
	profile(DDCENERGY, START);
   sys->energyInfo.called =0; 
	if (sys->potential[0] ->itype == ZEROPOTENTIAL && sys->npotential==0 )   // specical case of no interaction potential. Only kinetic terms needed. 
	{
		zeroAll(sys);
		if (e_eval_flag) 
		{
			kinetic_terms(sys, 1);
			eval_energyInfo(sys);
		}
		return 0; 
	}
	ETYPE *energyInfo = &sys->energyInfo; 

	cutoffs(ddc,sys);
	//auxNeighbor_begin(sys->nion); //HACK
	int update=0; 
	if (e_eval_flag == -1) 
	{
		e_eval_flag =0; 
		update=-1;
	}

	// Could do this differently by saying something like update =
	// loadBalance->balanceFunction(sys, loadBalance->parms); and having the
	// balanceFunction return something meaningful.  This would allow the
	// balance function to offer an opinion about whether an update is
	// needed.  The problem is, I'm not sure what update is before we get
	// here, and I'm not sure I want to clobber it with something like a
	// zero just because the load balance decides no update is needed.
   //
   if ( TEST0(sys->loop, ddc->loadBalance->rate) )
   {
      ddc->loadBalance->balanceFunction(ddc->loadBalance, sys, ddc);
      update=1;
   }

   //makeImageParticles(sys,ddc->rInfluence); 
   profile(UPDATEALL, START);
   ddcUpdateAll(ddc, sys, energyInfo, update);    
   profile(UPDATEALL, END);
   zeroAll(sys);
   profile(BARRIER1, START);
   WAIT(0);
   profile(BARRIER1, END);
   profile(P_FORCE, START);
   system_pCalculate(sys);
   if (state->q != NULL && ddc->lastUpdate == sys->loop) { for (unsigned i = 0; i < sys->nion; i++) state->q[i] = ((ATOMTYPE_PARMS*) (state->species[i])->parm)->charge;}
   
   for (int i = 0; i < sys->npotential; i++) sys->potential[i]->eval_potential((void *)sys, (void *)sys->potential[i]->parms, &(sys->energyInfo));

   if ((sys->neighborTableType & NEIGHBORTABLE_FAT) > 0) fsumX(sys, &(sys->energyInfo.virial));
   profile(P_FORCE, END);
   if (getSize(0) > 1) 
   {
      profile(BARRIER2, START);
      idleThreadsCheck(sys, COMM_LOCAL);  
      if (parityCheckingBarrier(COMM_LOCAL) != 0) return -1;
      profile(BARRIER2, END);
      ddcUpdateForce(ddc,sys->pCalculate);  
      /*ddcUpdateAccum(ddc); */

      profile(BARRIER3, START);
      WAIT(0);
      profile(BARRIER3, END);
   }
   if (sys->energyInfo.pdiaV != 0.0) for (unsigned i=0;i<sys->nlocal;i++)  state->virial[i]+=sys->energyInfo.pdiaV/sys->nglobal; 
   if (e_eval_flag) 
   {
      kinetic_terms(sys, 1);
      eval_energyInfo(sys);
   }
   profile(DDCENERGY, END);

   return 0;
}
void ddckinetic(int mode, SYSTEM*system)
{
   double rk_local;
   ETYPE *e=&system->energyInfo; 
   profile(KINETIC_TERMS, START);
   kinetic_terms(system, mode);
   rk_local = e->rk;
   MPI_Allreduce(&rk_local, &e->rk, 1, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   e->temperature = 2.0*e->rk/(3.0*system->nglobal);
   profile(KINETIC_TERMS, END);
}



/* Flops = 22 (MPM) */ 
/*
   void fker(int i, int j, double r2, double x, double y, double z)
   {
   double r,dr, fxij, fyij, fzij,v2;
   return; 
   r = sqrt(r2); 
   v2 = pair_function(i,j,r,&dr);
   fxij = -dr*x;
   fyij = -dr*y;
   fzij = -dr*z;
   fx[i] += fxij;
   fy[i] += fyij;
   fz[i] += fzij;
   fx[j] -= fxij;
   fy[j] -= fyij;
   fz[j] -= fzij;
   energy[i] += 0.5*v2;
   energy[j] += 0.5*v2;
   virial[i] += 0.1666666666666667*(fxij*x+fyij*y+fzij*z);
   virial[j] += 0.1666666666666667*(fxij*x+fyij*y+fzij*z);
   sion->xx -= fxij*x;
   sion->xy -= fxij*y;
   sion->xz -= fxij*z;
   sion->yy -= fyij*y;
   sion->yz -= fyij*z;
   sion->zz -= fzij*z;
   }
   */

/* Flops = 18  (MPM)   */ 
void fsumX(SYSTEM*sys, THREE_SMATRIX*restrict virial)
{
   NPARTICLE* particles = sys->neighbor->particles;
   double* fx = sys->collection->state->fx;
   double* fy = sys->collection->state->fy;
   double* fz = sys->collection->state->fz;
   double* energy = sys->collection->state->potentialEnergy;
   double* stateVirial = sys->collection->state->virial;

   for (unsigned i = 0; i < sys->nion; i++)
   {
      PAIRS* pij = particles[i].ifirst[INDEX];
      //		egi = &sys->collection->state->group[i]->e; 
      while (pij != NULL)
      {
         RPAIR* p = pij->p;
         if (p != NULL) 
         {
            int j = pij->j;
            //				egj = &sys->collection->state->group[j]->e; 
            double fxij = p->fp.x;
            double fyij = p->fp.y;
            double fzij = p->fp.z;
            double frxx = fxij*p->x; 
            double fryy = fyij*p->y; 
            double frzz = fzij*p->z; 
            //  			double frxy = fxij*p->y; 
            // 			double frxz = fxij*p->z;  
            //   			double fryz = fyij*p->z;  
            double v2 = 0.5*p->e;
            fx[i] += fxij;
            fy[i] += fyij;
            fz[i] += fzij;
            fx[j] -= fxij;
            fy[j] -= fyij;
            fz[j] -= fzij;
            //printf("%i fxi %f\n",i, fx[i]);
            //printf("%i fxj %f\n",j, fx[j]);
            energy[i] += v2;
            energy[j] += v2;
            stateVirial[i] += 0.1666666666666667*(frxx+fryy+frzz);
            stateVirial[j] += 0.1666666666666667*(frxx+fryy+frzz);

            virial->xx += fxij*p->x;
            virial->xy += fxij*p->y;
            virial->xz += fxij*p->z;
            virial->yy += fyij*p->y;
            virial->yz += fyij*p->z;
            virial->zz += fzij*p->z;
            /*
               if (egi == egj ) 
               {
               egi->eion += 2.0*v2;
               egi->virial.xx += frxx;
               egi->virial.xy += frxy;
               egi->virial.xz += frxz;
               egi->virial.yy += fryy;
               egi->virial.yz += fryz;
               egi->virial.zz += frzz;
               }
               else 
               {

               egi->f.x += fxij;
               egi->f.y += fyij;
               egi->f.z += fzij;

               egj->f.x -= fxij;
               egj->f.y -= fyij;
               egj->f.z -= fzij;

               egi->eion += v2;
               egi->virial.xx += (0.5*frxx);
               egi->virial.xy += (0.5*frxy);
               egi->virial.xz += (0.5*frxz);
               egi->virial.yy += (0.5*fryy);
               egi->virial.yz += (0.5*fryz);
               egi->virial.zz += (0.5*frzz);

               egj->eion += v2;
               egj->virial.xx += (0.5*frxx);
               egj->virial.xy += (0.5*frxy);
               egj->virial.xz += (0.5*frxz);
               egj->virial.yy += (0.5*fryy);
               egj->virial.yz += (0.5*fryz);
               egj->virial.zz += (0.5*frzz);
               }
               */

         }
         pij = pij->ilink;
      }
   }
}

/** This function does *not* reset the parity status.  The idea is that
 * some function higher up the call stack will recheck the status,
 * notice the problem and take corrective action.  At this point we just
 * want to notice the problem so that we can start to move up the call
 * stack.  */
int parityCheckingBarrier(MPI_Comm comm)
{
#ifndef BGL
   MPI_Barrier(comm);
   return 0;
#else
   profile(PARITYCHECK, START);
   rts_dcache_evict_transient();
   rts_dcache_evict_normal();
   int globalErr = 0;
   int localErr = parityStatus();
   MPI_Allreduce(&localErr, &globalErr, 1, MPI_INT, MPI_SUM, comm);
   profile(PARITYCHECK, END);
   return globalErr;
#endif	

}

void idleThreadsCheck(SYSTEM*system, MPI_Comm comm)
{
   int idleNodesSum = 0;
#ifdef WITH_OMP
   int idleNode = (system->nIdleThreads > 0);
   MPI_Allreduce(&idleNode, &idleNodesSum, 1, MPI_INT, MPI_SUM, comm);
#endif
   if ((idleNodesSum > 0) && (getRank(0) == 0)) 
      printf("Warning: %10d nodes have idle threads.\n", idleNodesSum);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
