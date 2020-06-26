#include "energyInfo.h"
#include <math.h>
#include "system.h"
#include "object.h"
#include "ptiming.h"
#include "mpiUtils.h"
#include "ddcMalloc.h"
#include "external.h"
static void allreduce(ETYPE *energyInfo)
{
	double psum[24], sum[24];
   int n=0;
	psum[n++] = energyInfo->rk;
	psum[n++] = energyInfo->eion;
	psum[n++] = energyInfo->virial.xx;
	psum[n++] = energyInfo->virial.yy;
	psum[n++] = energyInfo->virial.zz;
	psum[n++] = energyInfo->virial.xy;
	psum[n++] = energyInfo->virial.xz;
	psum[n++] = energyInfo->virial.yz;
	psum[n++] = energyInfo->pdiaV;
	psum[n++]  = energyInfo->f.x;
	psum[n++] = energyInfo->f.y;
	psum[n++] = energyInfo->f.z;
	psum[n++] = energyInfo->number;
	psum[n++] = energyInfo->e0_at_V;
	psum[n++] = energyInfo->thermal_flux.x; 
	psum[n++] = energyInfo->thermal_flux.y; 
	psum[n++] = energyInfo->thermal_flux.z; 
	psum[n++] = energyInfo->eBath; 
	psum[n++] = energyInfo->tion.xx;
	psum[n++] = energyInfo->tion.yy;
	psum[n++] = energyInfo->tion.zz;
	psum[n++] = energyInfo->tion.xy;
	psum[n++] = energyInfo->tion.xz;
	psum[n++] = energyInfo->tion.yz;
	MPI_Allreduce(psum, sum, n, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   n=0;
	energyInfo->rk = sum[n++];
	energyInfo->eion = sum[n++];
	energyInfo->virial.xx = sum[n++];
	energyInfo->virial.yy = sum[n++];
	energyInfo->virial.zz = sum[n++];
	energyInfo->virial.xy = sum[n++];
	energyInfo->virial.xz = sum[n++];
	energyInfo->virial.yz = sum[n++];
	energyInfo->pdiaV = sum[n++];
	energyInfo->f.x = sum[n++];
	energyInfo->f.y = sum[n++];
	energyInfo->f.z = sum[n++];
	energyInfo->number = sum[n++];
	energyInfo->e0_at_V = sum[n++];
	energyInfo->thermal_flux.x = sum[n++];
	energyInfo->thermal_flux.y = sum[n++];
	energyInfo->thermal_flux.z = sum[n++];
	energyInfo->eBath = sum[n++];
	energyInfo->tion.xx = sum[n++];
	energyInfo->tion.yy = sum[n++];
	energyInfo->tion.zz = sum[n++];
	energyInfo->tion.xy = sum[n++];
	energyInfo->tion.xz = sum[n++];
	energyInfo->tion.yz = sum[n++];
}
static unsigned testRates(SIGNED64 loop, unsigned *rates)
{
   for (unsigned i=1;i<=rates[0];i++)
   {
      if (loop%rates[i] == 0) 
      {
         return 1; 
      }
   }
   return 0; 
}
void eval_energyInfo(SYSTEM *sys)
{
   double vol; 
   ETYPE *eg; 
   profile(EVAL_ETYPE, START);
   ETYPE *energyInfo = &sys->energyInfo; 
   double vol_per_atom = 0.0; 
   if (testRates(sys->loop,sys->energyInfoObj->globalEvalRate))
   {
   vol = sys->box->volume;
   vol_per_atom = energyInfo->vol_per_atom; 
   if (vol_per_atom > 0.0) vol = vol_per_atom*energyInfo->number;
   else vol_per_atom = vol/energyInfo->number; 

   energyInfo->eion_local=energyInfo->eion;
   energyInfo->rk_local=energyInfo->rk;
   THREE_SMATRIX virial_local= energyInfo->virial;
   THREE_SMATRIX tion_local = energyInfo->tion;

   energyInfo->sion_local=virial_local;
   SMATACUM(energyInfo->sion_local, tion_local); 
   SMATNORM(energyInfo->sion_local, -vol/getSize(0));
   energyInfo->pion_local=-TRACE(energyInfo->sion_local)/3.0; 

   energyInfo->eBath=0.0; 
   for (int i=0;i<sys->ngroup;i++) 
   {
      GROUP *group=sys->group[i]; 
      group->energyBath=0.0; 
      if (group->energyBathFunction != NULL) group->energyBath=group->energyBathFunction(group); 
      energyInfo->eBath+=group->energyBath; 
   }
   allreduce(energyInfo);
   energyInfo->sion = energyInfo->virial; 
   SMATACUM(energyInfo->sion, energyInfo->tion); 
   SMATNORM(energyInfo->sion, -vol);

   VSCALE(energyInfo->thermal_flux,1.0/vol); 
   //THREE_SMATRIX sion=energyInfo->sion; 
   energyInfo->pion=-TRACE(energyInfo->sion)/3.0; 
   energyInfo->temperature = 2.0*energyInfo->rk/(3.0*energyInfo->number-sys->nConstraints);
   sys->energy = energyInfo->eion+energyInfo->rk; 
   }
   if (testRates(sys->loop,sys->energyInfoObj->groupEvalRate))
   {
      for (int i=0;i<sys->ngroup;i++) 
      {
         GROUP *group=sys->group[i]; 
         eg = &(group->energyInfo);
         eg->eBath=group->energyBath;
         allreduce(eg); 
         if (eg->number > 0.0)
         {
            vol = vol_per_atom*eg->number; 
            if (eg->vol_per_atom > 0.0) vol = eg->vol_per_atom*eg->number;
            eg->eion +=  energyInfo->e0_at_V*eg->number/energyInfo->number;
            //THREE_SMATRIX virial= eg->virial;
            //THREE_SMATRIX tion = eg->tion;
            eg->sion = eg->virial; 
            SMATACUM(eg->sion, eg->tion); 

            eg->sion.xx += energyInfo->pdiaV*eg->number/energyInfo->number;
            eg->sion.yy += energyInfo->pdiaV*eg->number/energyInfo->number;
            eg->sion.zz += energyInfo->pdiaV*eg->number/energyInfo->number;
            SMATNORM((eg->sion), -vol);
            eg->pion=-TRACE(eg->sion)/3.0; 
            eg->temperature = 2.0*eg->rk/(3.0*eg->number);
         }
      }	
   }

   profile(EVAL_ETYPE, END);
   profile(NOFIELD, AVE);
}
ENERGYINFO *energyInfo_init(void *parent, char *name)
{
   ENERGYINFO *energyInfoObj = (ENERGYINFO *) object_initialize(name, "ENERGYINFO", sizeof(ENERGYINFO));
   energyInfoObj->parent = parent; 
   energyInfoObj->groupEvalRate = ddcMalloc(17*sizeof(int)); 
   energyInfoObj->globalEvalRate = ddcMalloc(17*sizeof(int)); 
   energyInfoObj->groupEvalRate[0] = object_get((OBJECT *) energyInfoObj, "groupEvalRate",    energyInfoObj->groupEvalRate+1,   INT, 16, "1");
   energyInfoObj->globalEvalRate[0] = object_get((OBJECT *) energyInfoObj, "globalEvalRate",    energyInfoObj->globalEvalRate+1,   INT, 16, "1");
   return energyInfoObj; 

}
