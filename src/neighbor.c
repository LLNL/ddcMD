#include "neighbor.h"

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include "three_algebra.h"
#include "object.h"
#include "ptiming.h"
#include "error.h"
#include "system.h"
#include "ddcMalloc.h"
#include "heap.h"
#include "expandbuffer.h"
#include "preduce.h"
#include "units.h"

#define GI(i)  (*(gid_type *)(((char *)particleset->global_index)+i*particleset->stridelong64))
/*#define INDEX 1  */  /* PAIR FUNCTION OPT */
#define INDEX 0 

void pairlist(NBR *nbr, PARTICLESET*particleset);
void pairlist1(NBR *nbr, PARTICLESET*particleset);
void pairlist2(NBR *nbr, PARTICLESET*particleset);
void box_get(BOX_STRUCT *box ,int cmd, void *ptr);

void fker(int i,int j,double r1,double x,double y,double z);
int getRank(int); 

static NBR *nbr_static;
//static int nlocal = -1; 
extern	FILE *ddcfile; 

NBR *neighbor_init(void *parent, char *name)
{
	char *type; 
	NBR *nbr; 
	GEOM *geom; 
    nbr_static=nbr=(NBR *)object_initialize(name,"NEIGHBOR",sizeof(NBR));
	nbr->parent = parent; 
	nbr->particles = NULL; 
	nbr->nlocal=-1; 
	nbr->number_particles=-1; 
	nbr->npairs=-1; 
	nbr->rcut = NULL;
	nbr->nc = 0;
	nbr->Rpairs = NULL;
	nbr->Rpairs_blk = 0;
	object_get ((OBJECT*)nbr,"type",  &(type),STRING,1,"NORMAL");
	object_get((OBJECT *)nbr,"deltaR",&nbr->deltaR,WITH_UNITS,1,"0","l",NULL);
	nbr->geom=geom=geom_init(0.0, box_getBox(NULL)); // Note rcut argument is zero. Correct value of geom->rcut is set by neighborCutoff.
   nbr->lastUpdate= NO_NBR_UPDATE;

	object_get((OBJECT *)nbr,"minBoxSide",&geom->minBoxSide,WITH_UNITS,1,"0","l",NULL);
	return nbr; 
}
void neighborCutoff(NBR *nbr,int n, RCUT_TYPE *rcut)
{
   enum  type { small, new, exist}; 
	int k;
	if (n == -1)
	{
		nbr->nc =0; 
		return ; 
	}
   if (n == 0) return; 
	nbr->rcut = (RCUT_TYPE *) ExpandBuffers((void *)nbr->rcut, sizeof(RCUT_TYPE), nbr->nc+n, 16, LOCATION("neighborCutoff"),"nbr->rcut");
	for (int j=0;j<n;j++)
	{
      int type = small; 
      int i; 
		for (i=0;i<nbr->nc;i++)
		{
			if (nbr->rcut[i].value <rcut[j].value)  type = new ; 
			if (nbr->rcut[i].value==rcut[j].value)  type = exist; 
         if (type != small) break; 
		}
      if (type == small)  nbr->rcut[nbr->nc] = rcut[j]; 
      if (type ==new )
      {
		   for (k=nbr->nc;k>i;k--) nbr->rcut[k]=nbr->rcut[k-1];
		   nbr->rcut[i] = rcut[j]; 
      }
      if (type == exist)
      {
      }
		if (type != exist) nbr->nc++; 
	}
	int mode =0; 
	for (int i=0;i<nbr->nc;i++)
	{
		mode |= nbr->rcut[i].mode  ;
		nbr->rcut[i].mode = mode; 
	}

	nbr->geom->rcut = nbr->rcut[0].value + nbr->deltaR;
}

/** Determines whether the communication tables used to send remote
 *  particle information need to be updated.
 *
 *  The following conditions trigger an update:
 *
 *  0.  If the neighbor list hasn't yet been formed, i.e.,
 *  nbr->particles == NULL
 *
 *  1.  The number of local atoms obtained from the PATRICLESET has
 *  changed since the last call of this function.
 *
 *  2.  If the number of particles in the PATRICLESET and in nbr are
 *  different.
 *
 *  Returns 1 and sets nbr->update=1 if an update is required.  If no
 *  update is needed returns zero and leaves nbr->update
 *  unchanged. (Note that the value of nbr->update is never used.)
 */
int neighborCheck(NBR *nbr, PARTICLESET*particleset)
{
	int stride;
	double deltaR, x, y, z, dmax;
	THREE_VECTOR c,r,u={1.0,1.0,1.0},rbar; 
	GEOM* geom = nbr->geom;
	/*nbr->update=0;  */
	if (nbr->particles == NULL ) { nbr->update =1; return 1; }
	unsigned nlocal = *particleset->number_local;
	unsigned number_particles = *particleset->number_particles;
	if (nlocal  != nbr->nlocal) {nbr->update=1;return 1;}
	if (number_particles  != nbr->number_particles) {nbr->update=1;return 1;}
	/*if (fabs(nbr->rcut0-geom->rcut[0].value)/nbr->rcut0 > 0.1 )  {nbr->update=1;return 1;} */
	double* xptr = particleset->xptr;
	double* yptr = particleset->yptr;
	double* zptr = particleset->zptr;
	stride = particleset->stridedouble;
	deltaR = nbr->deltaR+(nbr->rcut[0].value-nbr->rcut0);
	dmax =0.0; 
	r = matrix_vector(*geom->hinv,u);
 	r = matrix_vector(nbr->h0,r);
	r.x = fabs(1.0-r.x); 
	r.y = fabs(1.0-r.y); 
	r.z = fabs(1.0-r.z); 
	dmax = r.x ; 
	if (dmax < r.y )  dmax = r.y  ;
	if (dmax < r.z )  dmax = r.z  ;
	dmax *=  (nbr->rcut[0].value+nbr->deltaR); 
	c =  nbr->center; 
	rbar.x = rbar.y = rbar.z = 0.0; 
	for (unsigned i = 0; i < nlocal; i++)
	{
		x = *(double *)(((char *)xptr) + i*stride)-c.x ;
		y = *(double *)(((char *)yptr) + i*stride)-c.y ;
		z = *(double *)(((char *)zptr) + i*stride)-c.z ;
		nearestImage(&x, &y, &z);  
		rbar.x += x+c.x ;
		rbar.y += y+c.y ;
		rbar.z += z+c.z ;
	}
	rbar.x /= nlocal;
	rbar.y /= nlocal;
	rbar.z /= nlocal;
	double d2max = 0.0;
   //int imax = -1; 
	for (unsigned i = 0; i < number_particles; i++)
//	for (unsigned i = 0; i < nlocal; i++)
	{
		THREE_VECTOR r1,r0;
		double dr2; 
		r1.x = ((*(double *)(((char *)xptr) + i*stride)))-rbar.x;
		r1.y = ((*(double *)(((char *)yptr) + i*stride)))-rbar.y;
		r1.z = ((*(double *)(((char *)zptr) + i*stride)))-rbar.z;
		nearestImage_fast(&r1.x, &r1.y, &r1.z);

		r0.x=nbr->particles[i].r0.x - nbr->rbar.x;
		r0.y=nbr->particles[i].r0.y - nbr->rbar.y;
		r0.z=nbr->particles[i].r0.z - nbr->rbar.z;
		nearestImage_fast(&r0.x, &r0.y, &r0.z);
 		//THREE_VECTOR u0 = matrix_vector(nbr->h0_inv,r0);
 		//THREE_VECTOR u1 = matrix_vector(*geom->hinv,r1);
		//x = u1.x - u0.x; if ( x > 0.5 ) x -= 1.0;  if ( x < -0.5) x += 1.0; 
		//y = u1.y - u0.y; if ( y > 0.5 ) y -= 1.0;  if ( y < -0.5) y += 1.0;
		//z = u1.z - u0.z; if ( z > 0.5 ) z -= 1.0;  if ( z < -0.5) z += 1.0;
      //
      x = r1.x-r0.x; 
      y = r1.y-r0.y; 
      z = r1.z-r0.z; 
		nearestImage_fast(&x, &y, &z);

		dr2 = (x*x+y*y+z*z); 
      if (dr2 > d2max)
      {
            d2max = dr2; 
      //    imax = i; 
      }
	}
	//CTHREE_VECTOR ev = eigenvalues(nbr->h0,&status);
	//d2max = d2max * (ev.x *ev.x) ; 
	dmax += 2.0*sqrt(d2max); 
	if (dmax < deltaR) return 0;
   //THREE_VECTOR rr; 
	//rr.x = ((*(double *)(((char *)xptr) + imax*stride)))-rbar.x;
	//rr.y = ((*(double *)(((char *)yptr) + imax*stride)))-rbar.y;
	//rr.z = ((*(double *)(((char *)zptr) + imax*stride)))-rbar.z;
   //THREE_VECTOR r0 = nbr->rbar;
   //THREE_VECTOR r1 = rbar;
   //double lc = units_convert(1.0, NULL, "Angstrom");
   //printf("%d: label= %llu r= %f %f %f r0= %f %f %f r1= %f %f %f dmax=%f deltaR=%f\n",getRank(0),GI(imax),lc*rr.x,lc*rr.y,lc*rr.z,lc*r0.x,lc*r0.y,lc*r0.z,lc*r1.x,lc*r1.y,lc*r1.z,lc*dmax,lc*deltaR); fflush(stdout); 
	nbr->update=1; 
	return 1; 
}
void neighborRef(NBR *nbr, PARTICLESET*particleset)
{
	double *xptr = particleset->xptr;
	double *yptr = particleset->yptr;
	double *zptr = particleset->zptr;
	int stride = particleset->stridedouble;
	nbr->number_particles = *particleset->number_particles;
	unsigned nlocal=nbr->nlocal = *particleset->number_local;
	nbr->particles = (NPARTICLE *) ExpandBuffers((void *)nbr->particles, sizeof(NPARTICLE), nbr->number_particles, 1024, LOCATION("neighborsRef"),"nbr->particles");
	nbr->rcut0 = nbr->rcut[0].value;
	
	THREE_VECTOR c =  nbr->center; 
	THREE_VECTOR rbar = {0.0,0.0,0.0};  
	for (unsigned  i = 0; i < nbr->number_particles; i++)
	{
		double x = *(double *)(((char *)xptr) + i*stride)-c.x ;
		double y = *(double *)(((char *)yptr) + i*stride)-c.y ;
		double z = *(double *)(((char *)zptr) + i*stride)-c.z ;
		nearestImage(&x, &y, &z);  
		nbr->particles[i].r0.x = x+c.x;
		nbr->particles[i].r0.y = y+c.y;
		nbr->particles[i].r0.z = z+c.z;
      if (i < nlocal) 
      {
		rbar.x += x+c.x ;
		rbar.y += y+c.y ;
		rbar.z += z+c.z ;
      }
	}
	rbar.x /= nlocal; 
	rbar.y /= nlocal;
	rbar.z /= nlocal;
	GEOM* geom = nbr->geom;
	nbr->h0= *geom->h; 
	nbr->h0_inv= *geom->hinv; 
	nbr->rbar = rbar; 
	return; 
}

/** This function updates the neighbor list nbr according to the data in
 *  particleset.
 *
 *  The PARTICLESET structure is mainly a convenient grouping of
 *  pointers to the actual atomic position data.
 *
 *  It isn't at all clear why this routine returns nbr.
 */
NBR *neighbors(NBR *nbr, PARTICLESET*particleset)
{
	profile(VBOX, START);
	GeomBox(nbr->geom, particleset);
	profile(VBOX, END);
	profile(PAIRLIST, START);
	pairlist(nbr, particleset);
	profile(PAIRLIST, END);
	return nbr;
}
NBR *neighbors1(NBR *nbr, PARTICLESET*particleset)
{
	profile(VBOX, START);
	GeomBox(nbr->geom, particleset);
	profile(VBOX, END);
	profile(PAIRLIST, START);
	if (nbr->geom->method ==  GEOM_P7N1ORTHO ) pairlist2(nbr, particleset);
	else pairlist1(nbr, particleset);
	profile(PAIRLIST, END);
	return nbr;
}
NBR *neighbors_eam_opt(NBR *nbr, PARTICLESET*particleset)
{
	void eam_opt_transform(void *);
	void eam_opt_box(void *);
	void eam_opt_neighbors(void *);
	return nbr; 
	eam_opt_transform(nbr->parent);
	profile(VBOX, START);
	eam_opt_box(nbr->parent);
	profile(VBOX, END);
	profile(PAIRLIST, START);
	eam_opt_neighbors(nbr->parent);
	profile(PAIRLIST, END);
	return nbr;
}
void pairelement(NBR* nbr, NPARTICLE *pi, PAIRS *pij, RPAIR *p, double x, double y, double z,double r2,int local_remote_flag)
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
	}
}
void pairUpdate(NBR *nbr, PARTICLESET*particleset)
{
	double* xptr = particleset->xptr;
	double* yptr = particleset->yptr;
	double* zptr = particleset->zptr;
	int stride = particleset->stridedouble;
	RCUT_TYPE* rcut = nbr->rcut;
	double rmax2 = rcut[0].value*rcut[0].value;
	nbr->update=0; 
	double minspan;
	box_get(NULL,MINSPAN,&minspan); 
	double R2cut = 0.25*minspan*minspan;
	int nRpair =0; 
	if (nbr->Rpairs != NULL) heapFree(nbr->Rpairs_blk);
	nbr->Rpairs=heapGet(&nbr->Rpairs_blk); 
	RPAIR* Rpairs = nbr->Rpairs; 
	PAIRS *pij=nbr->pairs+nbr->npairs; 
	int local_remote_flag =2; 
	
	int  nlocal=nbr->nlocal ;
	for (int i = nbr->number_particles-1; i>=0; i--)
	{
		NPARTICLE *pi = nbr->particles + i;
		THREE_VECTOR r; 
		
		if ( i+1 == nlocal  ) local_remote_flag=1; 
		for (int l=0;l<nbr->nc;l++)  pi->ifirst[l] =pi->ilast[l] = NULL; 
		r.x = *(double *)(((char *)xptr) + i*stride);
		r.y = *(double *)(((char *)yptr) + i*stride);
		r.z = *(double *)(((char *)zptr) + i*stride);
		for (int k = 0;k<pi->nb;k++)
		{
			pij--; 
			int j=pij->j ;
			double x = r.x - *(double *)(((char *)xptr) + j*stride) ;
			double y = r.y - *(double *)(((char *)yptr) + j*stride) ;
			double z = r.z - *(double *)(((char *)zptr) + j*stride) ;
			double r2 = x*x + y*y + z*z;
			if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
			local_remote_flag=3;  
			if (r2 < rmax2)  
			{
				pairelement(nbr,pi,pij,Rpairs+nRpair,x,y,z,r2,local_remote_flag);
				nRpair++; 
			}
			else
			{
				pij->p = NULL; 
				pij->ilink=NULL; 
			}
		}
	}
	for (unsigned i = 0; i < nbr->number_particles;  i++)
	{
		NPARTICLE *pi = nbr->particles + i;
		for (int l=nbr->nc-1;l>0;l--)
		{
			if ((pij = pi->ilast[l-1]) != NULL) pij->ilink = pi->ifirst[l]; 
			if (pi->ifirst[l-1] == NULL ) pi->ifirst[l-1] = pi->ifirst[l]; 
		}
	}
	nbr->nRpair=nRpair; 
	heapEndBlock(nbr->Rpairs_blk, nRpair*sizeof(RPAIR));
}
void pairelement1(int nc, double *rcut, NPARTICLE *pi, PAIRS *pij, RPAIR *pNew,  RPAIR *pOld)
{
	int l;
	for (l=1;l<nc;l++) {if (pOld->r > rcut[l]) break; }
	l--; 
	pij->ilink = pi->ifirst[l];
	pi->ifirst[l] = pij;
	if (pi->ilast[l] == NULL) pi->ilast[l] = pij ; 
}
void print_message(char *string) 
{
	if (getRank(0)==7)
	{
		int i; 
		PAIRS *pij; 
		NPARTICLE *pi; 
		i=77;
		pi = nbr_static->particles + i;
		/*if (pi->ifirst[0]->j > 1024)*/
		{
		pij = pi->ifirst[0];printf("%d: %d 0x%8p %d ",getRank(0),i,pij,pij->j); 
		pij = pi->ifirst[1]; printf("0x%8p %d ",pij,pij->j); 
		pij = pi->ifirst[2]; printf("0x%8p %d %s\n",pij,pij->j,string); 
		}
	}
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
