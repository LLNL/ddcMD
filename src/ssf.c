#include <stdio.h>
#include <complex.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

#include "analysis.h"
#include "realsizec.h"
#include "ddcMalloc.h"
#include "simulate.h"
#include "object.h"
#include "error.h"
#include "expandbuffer.h"

typedef struct ssf_parms_st
{
   unsigned   nkvectors;
	THREE_INT *n; 
	double kmax; 
} SSF_PARMS;

int getRank(int);
THREE_VECTOR kDirectionMax(THREE_VECTOR kbasis[3]);
static void  ssf_kvector(SSF_PARMS *parms);


SSF_PARMS* ssf_parms(ANALYSIS *analysis)
{
   SSF_PARMS* parms = (SSF_PARMS*) ddcMalloc(sizeof(SSF_PARMS)); 
	
   OBJECT* obj = (OBJECT*) analysis;

   object_get(obj, "kmax", &parms->kmax,DOUBLE,1,"1.0");
   
   parms->nkvectors = 0;
   parms->n = NULL;
   
   return parms;
}

void ssf_output(ANALYSIS* analysis)
{
}


void ssf_eval(ANALYSIS* analysis)
{
   SSF_PARMS* parms = (SSF_PARMS*) analysis->parms;
   ssf_kvector(parms);
   
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system;
   unsigned nlocal = sys->nlocal;
   double* rx = sys->collection->state->rx;
   double* ry = sys->collection->state->ry;
   double* rz = sys->collection->state->rz;
   

   THREE_VECTOR kbasis[3]; 
   box_get(sys->box, RECIP_LATTICEVECTORS, kbasis); 
   THREE_VECTOR Nmax = kDirectionMax(kbasis);

   // first build up the "basis"-vectors
   static double complex *px=NULL,*py=NULL,*pz=NULL;
	int nx_max  = 1 + parms->kmax * Nmax.x; 
	int ny_max  = 1 + parms->kmax * Nmax.y; 
	int nz_max  = 1 + parms->kmax * Nmax.z; 
	px = (double complex *)ddcRealloc(px, nlocal*nx_max*sizeof(*px)); 
	py = (double complex *)ddcRealloc(py, nlocal*ny_max*sizeof(*py));
	pz = (double complex *)ddcRealloc(pz, nlocal*nz_max*sizeof(*pz));
	
	for (unsigned i=0; i<nlocal; ++i)
	{
		THREE_VECTOR rv; 
		double complex* pxi = px+i; 
		double complex* pyi = py+i; 
		double complex* pzi = pz+i; 
		rv.x=rx[i]; rv.y=ry[i]; rv.z=rz[i]; 
		double complex cx = cexp(I*DOT(kbasis[0], rv));
		double complex cy = cexp(I*DOT(kbasis[1], rv));
		double complex cz = cexp(I*DOT(kbasis[2], rv));
		pxi[0] =  pyi[0] = pzi[0] =  1.0;  
		double complex sx,sy,sz; 
		sx = sy = sz = 1.0; 
		for (int m=1; m<nx_max; m++) pxi[m*nlocal] = (sx *= cx);
		for (int m=1; m<ny_max; m++) pyi[m*nlocal] = (sy *= cy);
		for (int m=1; m<nz_max; m++) pzi[m*nlocal] = (sz *= cz);
	}
	double fact[4]={1.0,0.5,0.25,0.0};
	for (unsigned l=0;l<parms->nkvectors;l++)
	{
		double complex ppp,ppm,pmp,pmm;
		double complex *pxi,*pyi,*pzi;
		double complex p[4],sp[4];
		unsigned nx,ny,nz;
		nx = parms->n[l].x; 
		ny = parms->n[l].y; 
		nz = parms->n[l].z; 
		double kx = nx*kbasis[0].x + ny*kbasis[1].x + nz*kbasis[2].x; 
		double ky = nx*kbasis[0].y + ny*kbasis[1].y + nz*kbasis[2].y; 
		double kz = nx*kbasis[0].z + ny*kbasis[1].z + nz*kbasis[2].z; 
		double k2 = kx*kx+ky*ky+kz*kz;
		pxi = px+(nlocal*nx);
		pyi = py+(nlocal*ny);
		pzi = pz+(nlocal*nz);
		ppp = ppm = pmp = pmm= 0.0;
		for (unsigned i=0;i<sys->nlocal;i++)
		{
			double complex qppp = pxi[i]*     pyi[i]*      pzi[i];
			double complex qppm = pxi[i]*     pyi[i] *conj(pzi[i]);
			double complex qpmp = pxi[i]*conj(pyi[i])*     pzi[i];
			double complex qpmm = pxi[i]*conj(pyi[i])*conj(pzi[i]);
			ppp += qppp;
			ppm += qppm;
			pmp += qpmp;
			pmm += qpmm;
		}
		p[0]=ppp;p[1]=ppm;p[2]=pmp;p[3]=pmm;
		MPI_Allreduce(p, sp, 8, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
		unsigned cnt_zero = (nx == 0) + (ny == 0) + (nz == 0) ;
		double Sk = fact[cnt_zero]*(sp[0]*conj(sp[0])+sp[1]*conj(sp[1])+sp[2]*conj(sp[2])+sp[3]*conj(sp[3]));
		printf("%f %f %f %f %e\n",k2,kx,ky,kz,Sk);
	}
	exit(0);

   
}

void ssf_close(ANALYSIS* analysis)
{
}

THREE_VECTOR kDirectionMax(THREE_VECTOR kbasis[3])
{	
	struct { double bb,bc,cb,cc;} Mi; 
	THREE_VECTOR V ,w; 
	double v[3]; 
	for (int i=0;i<3;i++)
	{
		THREE_VECTOR A=kbasis[(i+0)%3];
		THREE_VECTOR B=kbasis[(i+1)%3];
		THREE_VECTOR C=kbasis[(i+2)%3];
		double BB = DOT(B,B); 
		double CC = DOT(C,C); 
		double AB = DOT(A,B); 
		double AC = DOT(A,C); 
		double BC = DOT(B,C); 
		double det = BB*CC - BC*BC; 
		Mi.bb = CC/det; 
		Mi.bc = Mi.cb = -BC/det; 
		Mi.cc = BB/det; 
		double b = Mi.bb*AB + Mi.bc*AC;
		double c = Mi.cb*AB + Mi.cc*AC;
		V.x = A.x + b*B.x + c*C.x; 
		V.y = A.y + b*B.y + c*C.y; 
		V.z = A.z + b*B.z + c*C.z; 
		printf("%e %e %e\n",A.x,A.y,A.z); 
		v[i] = 1.0/sqrt(DOT(V,V)); 
	}
	w.x = v[0]; 
	w.y = v[1]; 
	w.z = v[2]; 
	return w; 
}
void ssf_kvector(SSF_PARMS *parms)
{
	int l,nkvector_est,cnt_zero;
	double  kx,ky,kz,k2,kmax;
	THREE_VECTOR kbasis[3];
	double vol = box_get_volume(NULL);
	box_get(NULL,RECIP_LATTICEVECTORS,kbasis); 
	THREE_VECTOR Nmax = kDirectionMax(kbasis);
	double kvol= 8.0*M_PI*M_PI*M_PI/vol; 
	kmax = parms->kmax; 
	double nmax = kmax /cbrt(kvol); 
	
	nkvector_est  = 1.2*(M_PI/6.0 * nmax *nmax * nmax  + 0.375 * M_PI * nmax*nmax) ;
	parms->n = (THREE_INT *) ExpandBuffers((void *)parms->n, sizeof(THREE_INT), nkvector_est, 64, LOCATION("ssf_kvector"),"parms->n");
	l=0; 
	int nx_max  = kmax * Nmax.x; 
	int ny_max  = kmax * Nmax.y; 
	int nz_max  = kmax * Nmax.z; 
	for (int nx=0;nx<nx_max;nx++)
	for (int ny=0;ny<ny_max;ny++)
	for (int nz=0;nz<nz_max;nz++) 
	{
		kx = nx*kbasis[0].x + ny*kbasis[1].x + nz*kbasis[2].x; 
		ky = nx*kbasis[0].y + ny*kbasis[1].y + nz*kbasis[2].y; 
		kz = nx*kbasis[0].z + ny*kbasis[1].z + nz*kbasis[2].z; 
		k2 = kx*kx+ky*ky+kz*kz; 
		if (k2 >= kmax*kmax) break ; 
		cnt_zero = (nx == 0) + (ny == 0) + (nz == 0) ;
		if ( cnt_zero < 3) 
		{
			if (    k2< 0.5*kmax*kmax ) 
			{
				if (l > nkvector_est) 
				{
					nkvector_est *= MAX(1.05,(1.05*(nx_max+1))/(nx+1.0));
					parms->n = (THREE_INT *) ExpandBuffers((void *)parms->n, sizeof(THREE_INT), nkvector_est, 64, LOCATION("ssf_kvector"),"parms->n");
				}
				parms->n[l].x=nx; 
				parms->n[l].y=ny; 
				parms->n[l].z=nz; 
				l++;
			}
		}
	}
	parms->nkvectors=l; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
