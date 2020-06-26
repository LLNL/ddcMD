#include <math.h>
#include <assert.h>
#include "preduce.h"
#include "back_communicate.h"
#include "bisectionLoadBalance.h"
#include "ddcMalloc.h"
void bisectionCalcParticleDestinations(DDC *ddc)
{
   static domain_info_p dinf = NULL;

   PARTICLE* ddcParticle = ddc->particles;
   BISECTION_BALANCE_PARMS *bbp = (BISECTION_BALANCE_PARMS *) (ddc->loadBalance->parms);
   assert(bbp != NULL);

   int np,pid;
   MPI_Comm_size(COMM_LOCAL,&np);
   MPI_Comm_rank(COMM_LOCAL,&pid);

   if(0) if(pid == 0) printf("In %s(): reassgin = %d\n",__func__,bbp->reassign);

   int bisection_success = 0;
   while (! bisection_success)
   {

      if (bbp->reassign == 1 || dinf == NULL)
      {
         /* Use bisection to update particle destination information */
         THREE_VECTOR c0v = box_get_corner(NULL),c1v = box_get_urcorner(NULL);
         double c0[3] = {c0v.x , c0v.y , c0v.z };
         double c1[3] = {c1v.x , c1v.y , c1v.z };

         assert(box_getType(NULL) == ORTHORHOMBIC);		 

         if (dinf != NULL) dinf_delete(dinf);
         dinf = redistribute2(np,pid,COMM_LOCAL, ddc->number_particles,ddcParticle, c0,c1);

         if (ddc->OverLapID != NULL)
         {
            ddcFree(ddc->OverLapID);
            ddc->OverLapID = NULL;
         }

         /* Set domain neighbors list. */
         {
            nodedata *neilist;
            const int nnei = dinf_get_neilist(dinf,ddc->rcut,&neilist);

            ddc->OverLapID = ddcRealloc(ddc->OverLapID, sizeof(int)*nnei);
            for (int i = 0; i<nnei; i++) ddc->OverLapID[i] = neilist[i].rank;
            ddc->nOverLap = nnei;
         }

         /* Set radius and center of local domain */
         {
            const nodedata *self = dinf_get_self(dinf);
            const double *c0 = self->c0,*c1 = self->c1;
            THREE_VECTOR c;
            double r;
            c.x = 0.5*(c0[0] + c1[0]);
            c.y = 0.5*(c0[1] + c1[1]);
            c.z = 0.5*(c0[2] + c1[2]);
            r = 0.0;
            for (int j = 0; j<3; j++)
            {
               const double d = 0.5*(c1[j] - c0[j]);
               r += d*d;
            }
            r = sqrt(r);
            domainset_setLocalCenter(&ddc->domains,c);
            domainset_setLocalRadius(&ddc->domains,r);
         }
         bbp->reassign = 0;
         bisection_success = 1;
      }
      else
      {
         nodedata *neilist;
         const nodedata *self;
         int nnei;
         int outside_error = 0;

         assert(dinf != NULL);

         bisection_success = 1;

         self = dinf_get_self(dinf);
         nnei = dinf_get_neilist(dinf,ddc->rcut,&neilist);
         for (unsigned int i = 0; i<ddc->number_particles; i++)
         {
            double pos[3];
            pos[0] = ddcParticle[i].r.x;
            pos[1] = ddcParticle[i].r.y;
            pos[2] = ddcParticle[i].r.z;
            backInBox(pos+0,pos+1,pos+2);
            if (dinf_inside(self->c0,self->c1,pos))
            {
               ddcParticle[i].domain_id = pid;
            }
            else
            {
               int j;
               for (j = 0; j<nnei; j++)
               {
                  if (dinf_inside(neilist[j].c0,neilist[j].c1,pos))
                  {
                     ddcParticle[i].domain_id = neilist[j].rank;
                     break;
                  }
               }
               if (j >= nnei)
               {
                  outside_error = 1;
                  break;
               }
            }
         }

         /* Rerun bisection from scratch if any task had particles outside its domain neighbor list */
         {
            int noutside;
            MPI_Allreduce(&outside_error,&noutside,1,MPI_INT,MPI_SUM,COMM_LOCAL);
            if (noutside > 0)
            {
               bbp->reassign = 1;
               bisection_success = 0;
               if (pid == 0)
               {
                  printf("In %s() at %s:%d: Warning: Rerunning bisection on non-loadbalance step, because %d "
                        "tasks had particles outside their domain neighbor lists.\n",
                        __func__,__FILE__,__LINE__,noutside);
               }
            }
         }
      }
   }
}

