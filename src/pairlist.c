#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "ptiming.h"
#include "error.h"
#include "utilities.h"
#include "system.h"
#include "ddcMalloc.h"
#include "expandbuffer.h"
#include "neighbor.h"
#include "heap.h"
#include "preduce.h"
#define MIN_NPAIR_PER_PARTICLE  16 

int getRank(int);

#define PCHECK(x) x 
#define GI(i)  (*(gid_type *)(((char *)particleset->global_index)+i*particleset->stridelong64))
/*#define INDEX 1  */  /* PAIR FUNCTION OPT */
#define INDEX 0 
static void pairelementStatic (NBR* nbr, NPARTICLE *pi, PAIRS *pij, RPAIR *p, double x, double y, double z,double r2,int local_remote_flag)
{
	int l;
	RCUT_TYPE *rcut = nbr->rcut; 
	for (l=1;l<nbr->nc;l++) {if (r2 > rcut[l].value*rcut[l].value) break; }
	l--; 
	pij->ilink = pi->ifirst[l];
	pi->ifirst[l] = pij;
	if (pi->ilast[l] == NULL) pi->ilast[l] = pij ; 
	/*if ( (rcut[l].mode & local_remote_flag) != 0  )    */
	{
		pij->p = p ;//= Rpairs+nRpair;
		p->x = x;
		p->y = y;
		p->z = z;
		p->r = sqrt(r2); 
      //PCHECK(printf("init %i  %f\n",pij->j, p->fp.x);)
	}
}

void pairlist(NBR *nbr, PARTICLESET*particleset)
{
	double volume;
	int j, k, l, local_remote_flag;
	static int npairlast=0; 
	static int mpair = 0;
	static int npair_per_particle=MIN_NPAIR_PER_PARTICLE;
	static PAIRS *pairs = NULL;
	int nbi; 
	GEOMBOX *box_neigh, *box;
	THREE_VECTOR r;
	double minspan, r2, x, y, z;

	double* xptr = particleset->xptr;
	double* yptr = particleset->yptr;
	double* zptr = particleset->zptr;
	int stride = particleset->stridedouble;
	GEOM* geom = nbr->geom; 
	RCUT_TYPE* rcut = nbr->rcut;
	double deltaR=(nbr->deltaR);
	double rmax_plus_delta2 = ((rcut[0].value+deltaR)*(rcut[0].value+deltaR));
	double rmax2 = rcut[0].value*rcut[0].value;
	box_get(NULL,MINSPAN,&minspan); 
	double R2cut = 0.25*minspan*minspan;
	box_get(NULL,VOLUME,&volume); 
	unsigned nion = *particleset->number_particles;
	unsigned nlocal = *particleset->number_local;
	npairlast = nbr->npairs; 
   int nSearch =0; 
	int npair = 0;
	int nRpair =0; 
	
	if (nbr->Rpairs != NULL)
		heapFree(nbr->Rpairs_blk);
	nbr->Rpairs=heapGet(&nbr->Rpairs_blk); 
	RPAIR* Rpairs=nbr->Rpairs; 
   //printf("rpair %f\n",Rpairs[0].fp.x);
	/*mpair = ExpandBuffersMaxLength(pairs)/sizeof(PAIRS); */
	if (mpair < npair_per_particle*nion * 1.10  ) 
	{
		mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
		ExpandBuffersMode("FREE_MALLOC");
		pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
		ExpandBuffersMode("REALLOC");
		assert(pairs != NULL) ; 
	}
	heapTestAvailable(nbr->Rpairs_blk,  sizeof(RPAIR)*mpair);
	local_remote_flag = 1;   
	unsigned n = nion; 
	for (unsigned i = 0; i < n; i++) nbr->particles[i].nb=0; 
	if (rcut[0].mode == RCUT_LOCAL ) n = nlocal;  
   THREE_VECTOR rv[nion]; 
   gid_type gi[nion]; 
	for (unsigned i=0;i<nion;i++) 
   { 
		rv[i].x = *(double *)(((char *)xptr) + i*stride) ;
		rv[i].y = *(double *)(((char *)yptr) + i*stride) ;
		rv[i].z = *(double *)(((char *)zptr) + i*stride) ;
      gi[i] = GI(i); 
      
   }
	for (unsigned i = 0; i < n; i++)
	{
		if (npair + 16*npair_per_particle > mpair) //If buffers not big enough start over. 
		{
			if(getRank(0)==0) printf("Expanding pairs: npair=%d mpair=%d\n",npair,mpair);
			npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
			mpair = MAX(mpair,npair + 32*npair_per_particle);
			mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
			heapTestAvailable(nbr->Rpairs_blk,  sizeof(RPAIR)*mpair);
			ExpandBuffersMode("FREE_MALLOC");
			pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
			assert(pairs != NULL) ; 
			ExpandBuffersMode("REALLOC");
			i=0; 
			nSearch = 0;
			npair = 0;
			nRpair =0; 
			Rpairs=nbr->Rpairs; 
			local_remote_flag = 1;  
		}
		NPARTICLE *pi = nbr->particles + i;
		for (l=0;l<nbr->nc;l++)  pi->ifirst[l] =pi->ilast[l] = NULL; 
		gid_type gi_i = gi[i];;
      r = rv[i]; 
		box = geom->pinfo[i].box;

		nbi = 0;
		if ( i == nlocal ) local_remote_flag = 2;   //This statement is ignored. Why?
      local_remote_flag=3;  
		PAIRS *pij = pairs+npair; 
		{
			for (k = 0; k < box->nn; k++)
			{
				box_neigh = box + box->nlist[k];
				j = box_neigh->first;
				while (j != -1)
				{
					if (gi_i < gi[j])
					{
						x = r.x - rv[j].x;
						y = r.y - rv[j].y;
						z = r.z - rv[j].z;
						//nearestImage_fast(&x, &y, &z);
                  nSearch++; 
						pij->p=NULL; 
						r2 = x*x + y*y + z*z;
						if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
						if (r2 < rmax_plus_delta2)
						{
							pij->j = j;
                     if (r2 < rmax2) 
                     {
                     	pairelementStatic(nbr,pi,pairs+npair,Rpairs+nRpair,x,y,z,r2,local_remote_flag);
/*
                        pij->ilink = pi->ifirst[0];
                        pi->ifirst[0] = pij;
                        if (pi->ilast[0] == NULL) pi->ilast[0] = pij ; 
                        RPAIR *p; 
                        pij->p = p = Rpairs+nRpair;
                        p->x = x;
                        p->y = y;
                        p->z = z;
                        p->r = sqrt(r2); 
*/
                        nRpair++;
                     }
                     npair++;
                     nbi++;
                     pij++; 
                  }
               }
               j = geom->pinfo[j].next;
            }
         }
      }
      pi->nb = nbi;
      npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
   }
   npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,npair/MAX(nion,1));
   for (k=npair;k<npairlast;k++) pairs[k].p=NULL;  
   for (unsigned i = 0; i < nbr->number_particles;  i++)
   {
      NPARTICLE *pi = nbr->particles + i;
      for (l=nbr->nc-1;l>0;l--)
      {
         PAIRS *pij; 
         if ((pij = pi->ilast[l-1]) != NULL ) pij->ilink = pi->ifirst[l]; 
         if (pi->ifirst[l-1] == NULL ) pi->ifirst[l-1] = pi->ifirst[l]; 
      }
   }
   nbr->update=0; 
   nbr->pairs = pairs; 
   nbr->nSearch=nSearch; 
   nbr->npairs=npair; 
   nbr->nRpair=nRpair; 
   heapEndBlock(nbr->Rpairs_blk, nRpair*sizeof(RPAIR));
   if (getRank(0)==-1) printf("npairs per particle=%f \n",(1.0*npair)/nion);

}
void pairlist1(NBR *nbr, PARTICLESET*particleset)
{
   NPARTICLE *pi;
   double volume;
   int j, k;
   gid_type gi_i, gi_j;
   static int mpair = 0;
   static int npair_per_particle=MIN_NPAIR_PER_PARTICLE;
   static PAIRS *pairs = NULL;
   GEOMBOX *box_neigh, *box;
   PAIRS *pij;
   THREE_VECTOR r;
   double minspan, r2, x, y, z;

   double* xptr = particleset->xptr;
   double* yptr = particleset->yptr;
   double* zptr = particleset->zptr;
   int stride = particleset->stridedouble;
   GEOM* geom = nbr->geom; 
   RCUT_TYPE* rcut = nbr->rcut;
   double deltaR=(nbr->deltaR);
   double rmax_plus_delta2 = ((rcut[0].value+deltaR)*(rcut[0].value+deltaR));
   box_get(NULL,MINSPAN,&minspan); 
   double R2cut = 0.25*minspan*minspan;
   box_get(NULL,VOLUME,&volume); 
   unsigned nion = *particleset->number_particles;
   unsigned nlocal = *particleset->number_local;
   int nSearch = 0;
   int npair = 0;
   if (mpair < npair_per_particle*nion * 1.10  ) 
   {
      mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
      ExpandBuffersMode("FREE_MALLOC");
      pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
      ExpandBuffersMode("REALLOC");
   }
   unsigned n = nion; 
   for (unsigned  i = 0; i < n; i++) nbr->particles[i].nb=0; 
   // calculate local-local pairs (single sided) local-remote and remote->local and remote-remote pairs. 
   if (rcut[0].mode == RCUT_LOCAL ) n = nlocal;   //  Skip calculating remote->local and remote-remote pairs. 
   for (unsigned i = 0; i < n; i++)
   {
      if (npair + 16*npair_per_particle > mpair) //Start over if buffers not large enough. 
      {
         if(getRank(0)==0) printf("Expanding pairs: npair=%d mpair=%d\n",npair,mpair);
         npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
         mpair = MAX(mpair,npair + 32*npair_per_particle);
         mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
         ExpandBuffersMode("FREE_MALLOC");
         pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
         if(getRank(0)==0) printf("new mpair=%d\n",mpair);
         ExpandBuffersMode("REALLOC");
         i=0; 
         nSearch = 0;
         npair = 0;
      }
      pi = nbr->particles + i;
      pi->ifirst[0] =pi->ilast[0] = NULL; 
      gi_i = GI(i);
      r.x = *(double *)(((char *)xptr) + i*stride) ;
      r.y = *(double *)(((char *)yptr) + i*stride) ;
      r.z = *(double *)(((char *)zptr) + i*stride) ;
      box = geom->pinfo[i].box;
      //start = pairs + npair;
      int nbi = 0;
      pij = pairs+npair; 
      {
         for (k = 0; k < box->nn; k++)
         {
            box_neigh = box + box->nlist[k];
            j = box_neigh->first;
            while (j != -1)
            {
               gi_j = GI(j);
               if (gi_i < gi_j)
               {
                  x = r.x - *(double *)(((char *)xptr) + j*stride) ;
                  y = r.y - *(double *)(((char *)yptr) + j*stride) ;
                  z = r.z - *(double *)(((char *)zptr) + j*stride) ;
                  nSearch++; 
                  pij->p=NULL; 
                  r2 = x*x + y*y + z*z;
                  if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
                  if (r2 < rmax_plus_delta2)
                  {
                     pij->j = j;
                     pij->ilink = pi->ifirst[0];
                     pi->ifirst[0] = pij;
                     if (pi->ilast[0] == NULL) pi->ilast[0] = pij ; 
                     npair++;
                     nbi++;
                     pij++; 
                  }
               }
               j = geom->pinfo[j].next;
            }
         }
      }
      pi->nb = nbi;
      npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
   }
   npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,npair/MAX(nion,1));
   nbr->update=0; 
   nbr->pairs = pairs; 
   nbr->nSearch=nSearch; 
   nbr->npairs=npair; 
   nbr->nRpair=0; 
   if (getRank(0)==-1) printf("npairs per particle=%f \n",(1.0*npair)/nion);

}
void pairlist2(NBR *nbr, PARTICLESET*particleset)
{
   NPARTICLE *pi;
   double volume;
   int j, k;
   gid_type gi_i, gi_j;
   static int mpair = 0;
   static int npair_per_particle=MIN_NPAIR_PER_PARTICLE;
   static PAIRS *pairs = NULL;
   GEOMBOX *box_neigh, *box;
   PAIRS *pij;
   THREE_VECTOR r;
   double minspan, r2, x, y, z;

   double* xptr = particleset->xptr;
   double* yptr = particleset->yptr;
   double* zptr = particleset->zptr;
   int stride = particleset->stridedouble;
   GEOM* geom = nbr->geom; 
   RCUT_TYPE* rcut = nbr->rcut;
   double deltaR=(nbr->deltaR);
   double rmax_plus_delta2 = ((rcut[0].value+deltaR)*(rcut[0].value+deltaR));
   box_get(NULL,MINSPAN,&minspan); 
   double R2cut = 0.25*minspan*minspan;
   box_get(NULL,VOLUME,&volume); 
   unsigned nion = *particleset->number_particles;
   unsigned nlocal = *particleset->number_local;
   int nSearch = 0;
   int npair = 0;
   if (mpair < npair_per_particle*nion * 1.10  ) 
   {
      mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
      ExpandBuffersMode("FREE_MALLOC");
      pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
      ExpandBuffersMode("REALLOC");
   }
   unsigned n = nion; 
   for (unsigned  i = 0; i < n; i++) nbr->particles[i].nb=0; 
   // calculate local-local pairs (single sided) local-remote and remote->local and remote-remote pairs. 
   if (rcut[0].mode == RCUT_LOCAL ) n = nlocal;   //  Skip calculating remote->local and remote-remote pairs. 
   for (unsigned i = 0; i < n; i++)
   {
      if (npair + 16*npair_per_particle > mpair) //Start over if buffers not large enough. 
      {
         if(getRank(0)==0) printf("Expanding pairs: npair=%d mpair=%d\n",npair,mpair);
         npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
         mpair = MAX(mpair,npair + 32*npair_per_particle);
         mpair = MAX(mpair,npair_per_particle*nion * 1.10 ); 
         ExpandBuffersMode("FREE_MALLOC");
         pairs = (PAIRS *) ExpandBuffers((void *)pairs, sizeof(PAIRS), mpair, 1024, LOCATION("pairlist"),"pairs");
         if(getRank(0)==0) printf("new mpair=%d\n",mpair);
         ExpandBuffersMode("REALLOC");
         i=0; 
         nSearch = 0;
         npair = 0;
      }
      pi = nbr->particles + i;
      pi->ifirst[0] =pi->ilast[0] = NULL; 
      gi_i = GI(i);
      r.x = *(double *)(((char *)xptr) + i*stride) ;
      r.y = *(double *)(((char *)yptr) + i*stride) ;
      r.z = *(double *)(((char *)zptr) + i*stride) ;
      box = geom->pinfo[i].box;
      //start = pairs + npair;
      int nbi = 0;
      pij = pairs+npair; 
      {
         for (k = 0; k < box->nn; k++)
         {
            box_neigh = box + box->nlist[k];
            int type= box_neigh->type; 
            j = box_neigh->first;
            if (type != 13)
            {
               int ix = box_neigh->ix; 
               int iy = box_neigh->iy; 
               int iz = box_neigh->iz; 
              //  printf("%2d %2d %2d %4d %2d\n",ix,iy,iz,j,type); 
               int nx = geom->nx; 
               int ny = geom->ny; 
               int nz = geom->nz; 

                ix = -(type%3-1); 
                iy = -((type/3)%3-1); 
                iz = -((type/9)%3-1); 
               // printf("%2d %2d %2d\n",ix,iy,iz); 
               
               int offset  = ix*(nx-2) + iy*nx*(ny-2) + iz*nx*ny*(nz-2); 
               //printf("offset = %d \n",offset); 
               box_neigh+=offset; 
               j = box_neigh->first;
             type= box_neigh->type; 
                ix = box_neigh->ix; 
                iy = box_neigh->iy; 
                iz = box_neigh->iz; 
                //printf("%2d %2d %2d %4d %2d\n\n",ix,iy,iz,j,type); 
            }
            while (j != -1)
            {
               gi_j = GI(j);
               if (gi_i < gi_j)
               {
                  x = r.x - *(double *)(((char *)xptr) + j*stride) ;
                  y = r.y - *(double *)(((char *)yptr) + j*stride) ;
                  z = r.z - *(double *)(((char *)zptr) + j*stride) ;
                  nSearch++; 
                  pij->p=NULL; 
                  r2 = x*x + y*y + z*z;
                  if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
                  if (r2 < rmax_plus_delta2)
                  {
                     pij->j = j;
                     pij->ilink = pi->ifirst[0];
                     pi->ifirst[0] = pij;
                     if (pi->ilast[0] == NULL) pi->ilast[0] = pij ; 
                     npair++;
                     nbi++;
                     pij++; 
                  }
               }
               j = geom->pinfo[j].next;
            }
         }
      }
      pi->nb = nbi;
      npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,1.0 + npair/(i+1.0));
   }
   printf("npair=%d\n",npair); 
   npair_per_particle = MAX(MIN_NPAIR_PER_PARTICLE,npair/MAX(nion,1));
   nbr->update=0; 
   nbr->pairs = pairs; 
   nbr->nSearch=nSearch; 
   nbr->npairs=npair; 
   nbr->nRpair=0; 
   if (getRank(0)==-1) printf("npairs per particle=%f \n",(1.0*npair)/nion);

}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
