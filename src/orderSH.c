#include "orderSH.h"
#include  <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include <unistd.h>
#include  <mpi.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "pio.h"
#include "crc32.h"
#include "simulate.h"
#include "system.h"
#include  "sph.h"
#include "ddcMalloc.h"
#include "expandbuffer.h"
#include "ioUtils.h"
#include "io.h"
#include "preduce.h"

#define LMAX  16
typedef struct complex_st
{
	double real, imag;
} COMPLEX;

typedef struct ylm_st
{
	COMPLEX value, x, y, z;
} YLM;

YLM q[LMAX];
typedef struct cluster_st
{
	unsigned label, glabel;
	unsigned int size;
	double Rrms;
	double radius;
	THREE_VECTOR Rave;
	COMPLEX qAve[LMAX];
} CLUSTER;

enum thresholds
{ NOGROUP = -1, LIQUID = 0, INTERFACE = 1, CRYSTAL = 2, HIGHORDER = 3 };
static double W;
static int L;
static COMPLEX Y[LMAX];
static COMPLEX **qlocal = NULL;
static COMPLEX *qAccum = NULL;
static double **qnorm = NULL;
static double *Q = NULL;
static double *Wlocal = NULL;
static int *C = NULL;
static int *G = NULL;
static int Lv[16]; 
static int nL;
static NPARTICLE *particles = NULL;
static int IndexR1 = -1;
static int IndexR2 = -1;
static MPI_Datatype CLUSTER_TYPE;
CLUSTER cluster[512];

FUNCTION *function_init(void *,char *);
void sph(SPH*p, double x, double y, double z, YLM*ylm);
SPH *sphinit(int);
int getRank(int);
int getSize(int);
size_t Pwrite(const void *ptr, size_t size, size_t nmemb, PFILE*file);

static double wfunc(ORDERSH_PARMS*p, double r, double *dwdr);
static void orderSHlocal(SYSTEM*sys, ORDERSH_PARMS*parms, ETYPE*e);
static double orderDot(COMPLEX*a, COMPLEX*b, int L);


ORDERSH_PARMS *orderSH_parms(OBJECT*object)
{
	ORDERSH_PARMS *parms;
	char *functionname;
	parms = ddcMalloc(sizeof(ORDERSH_PARMS));
	object_get((OBJECT *) object, "function", &functionname, STRING, 1, "LINEAR");
	parms->function = function_init(object,functionname);
	nL = object_get((OBJECT *) object, "L", Lv, INT, 16, "6");
	object_get((OBJECT *) object, "Vo", &parms->Vo, DOUBLE, 1, "0.0");
	parms->Lo = 0.0;
	if (parms->Vo == 0.0)
	{
	   object_get((OBJECT *) object, "Lo", &parms->Lo, WITH_UNITS, 1, "0.0","l",NULL);
	}
	object_get((OBJECT *) object, "r1o", &parms->r1o, WITH_UNITS, 1, "0.0","l",NULL);
	object_get((OBJECT *) object, "r2o", &parms->r2o, WITH_UNITS, 1, "0.0","l",NULL);
	object_get((OBJECT *) object, "lamda", &parms->lamda, WITH_UNITS, 1, "0.0","m*l^2/t^2",NULL);
	object_get((OBJECT *) object, "globalevalrate", &parms->globalevalrate, INT, 1, "-1");
	object_get((OBJECT *) object, "localevalrate", &parms->localevalrate, INT, 1, "-1");
	if (parms->localevalrate == -1 )  parms->localevalrate = ((SIMULATE *) (object->parent->parent))->snapshotrate;
	if (parms->globalevalrate == -1 ) parms->globalevalrate = ((SIMULATE *) (object->parent->parent))->printrate;
	parms->sph = ddcMalloc(nL*sizeof(SPH *));
	parms->nL = nL;
	qlocal = (COMPLEX **) ExpandBuffers((void *)qlocal, sizeof(COMPLEX *), parms->nL, 16, LOCATION("orderSH_parms"), "qlocal");
	qnorm = (double **)ExpandBuffers((void *)qnorm, sizeof(double *), parms->nL, 16, LOCATION("orderSH_parms"), "qnorm");
	for (int i = 0; i < nL; i++)
	{
		L = Lv[i];
		parms->sph[i] = sphinit(L);
		qlocal[i] = NULL;
		qnorm[i] = NULL;
	}
	parms->sph[0] = sphinit(Lv[0]);
	L = parms->L = Lv[0];
	{
		int i, n, blkcnt[10];
		MPI_Datatype types[10];
		MPI_Aint disp[10];
		n = 3;
		blkcnt[0] = 2;
		blkcnt[1] = 5;
		blkcnt[2] = 2*LMAX;
		types[0] = MPI_UNSIGNED;
		types[1] = MPI_DOUBLE;
		types[2] = MPI_DOUBLE;
		MPI_Get_address(&cluster[0].label, &disp[0]);
		MPI_Get_address(&cluster[0].Rrms, &disp[1]);
		MPI_Get_address(&cluster[0].qAve, &disp[2]);
		for (i = n; i >= 0; i--)
			disp[i] -= disp[0];
		MPI_Type_create_struct(n, blkcnt, disp, types, &CLUSTER_TYPE);
	}
	return parms;
}

RCUT_TYPE *orderSHCutoff(SYSTEM*sys, ORDERSH_PARMS*parms, int *n)
{
	static RCUT_TYPE rcut[4];
	static int ncut = 2;
	double vol, scale, l;
	scale = 1.0;
	if (parms->Vo > 0.0)
	{
		vol = sys->box->volume/sys->nglobal;
		scale = cbrt(vol/parms->Vo);
	}
	if (parms->Lo > 0.0)
	{
		vol = sys->box->volume;
		l = cbrt(vol);
		scale = l/parms->Lo;
	}
	ncut =2; 
	parms->r1 = parms->r1o*scale;
	parms->r2 = parms->r2o*scale;
	rcut[0].value = parms->r2;
	rcut[1].value = parms->r1;
	rcut[0].mode=RCUT_ALL; 
	rcut[1].mode=RCUT_ALL; 
	rcut[ncut] = rcut[0];
	*n = ncut;
	if (TEST0(sys->loop, parms->localevalrate)) rcut[ncut].value = 2.0*rcut[0].value;

	return rcut;
}

double wfunc(ORDERSH_PARMS*order, double r, double *dwdr)
{
	double w, dri;
	*dwdr = 0.0;
	if (r < order->r1) return 1.0;
	if (r > order->r2) return 0.0;
	dri = 1.0/(order->r2 - order->r1);
	w = 0.5 + 0.5*cos(M_PI*(r - order->r1)*dri);
	*dwdr = -0.5*dri*M_PI*sin(M_PI*(r - order->r1)*dri);
	return w;
}

double orderPass1(SYSTEM *sys, ORDERSH_PARMS*order)
{
	int i, m, nlocal;
	double phi, w, dwdr;
	COMPLEX psum[LMAX], sum[LMAX];
	double (*fcn) ();
	PAIRS *pij;
	RPAIR *p;
	nlocal = sys->nlocal;
	W = 0;
	L = order->sph[0]->L;
	for (m = 0; m <= L; m++) Y[m].real = Y[m].imag = 0.0;
	for (i = 0; i < nlocal; i++)
	{
		pij = particles[i].ifirst[IndexR2];
		while (pij  != NULL)
		{
			p =pij->p;
			w = wfunc(order, p->r, &dwdr);
			if (w > 0.0)
			{
				sph(order->sph[0], p->x, p->y, p->z, q);
				for (m = 0; m <= L; m++)
				{
					Y[m].real += q[m].value.real*w;
					Y[m].imag += q[m].value.imag*w;
				}
				W += w;
			}
			pij = pij->ilink;
		}
	}

	for (m = 0; m <= L; m++) psum[m] = Y[m];
	psum[L + 1].real = W;
	MPI_Allreduce(psum, sum, 2*L + 3, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
	for (m = 0; m <= L; m++) Y[m] = sum[m];
	W = sum[L + 1].real;

	Y[0].real /= W;
	phi = Y[0].real*Y[0].real;
	for (m = 1; m <= L; m++)
	{
		Y[m].real /= W;
		Y[m].imag /= W;
		phi += 2.0*(Y[m].real*Y[m].real + Y[m].imag*Y[m].imag);
	}
	phi *= 4.0*M_PI/(2.0*L + 1.0);
	fcn = order->function->f;
	order->phi = phi;
	order->Energy = fcn(order->function->parms, order->phi, &order->dEnergydphi);
	order->Energy *= order->lamda;
	order->dEnergydphi *= order->lamda;
	return nlocal*order->Energy;
}

void orderPass2(ORDERSH_PARMS*order, unsigned nlocal, unsigned nsimul)
{
	if (order->lamda == 0.0) return;
	L = order->sph[0]->L;
	double c = 4.0*nsimul*4.0*M_PI/(2.0*L + 1.0)*order->dEnergydphi/W;
	c=0; 
	if (c < 1e-16) return;
	for (unsigned i = 0; i < nlocal; i++)
	{
		PAIRS* pij = particles[i].ifirst[IndexR2];
		while (pij  != NULL)
		{
			RPAIR* p = pij->p;
			double dwdr;
			double w = wfunc(order, p->r, &dwdr);
			if (w > 0.0)
			{
				dwdr /= (p->r);
				double x = p->x;
				double y = p->y;
				double z = p->z;
				sph(order->sph[0], x, y, z, q);
				int m = 0;
				p->fp.x -= 0.5*c*Y[m].real*(q[m].x.real*w + x*dwdr*(q[m].value.real - Y[m].real));
				p->fp.y -= 0.5*c*Y[m].real*(q[m].y.real*w + y*dwdr*(q[m].value.real - Y[m].real));
				p->fp.z -= 0.5*c*Y[m].real*(q[m].z.real*w + z*dwdr*(q[m].value.real - Y[m].real));
				for (m = 1; m <= L; m++)
				{
					p->fp.x -= c*Y[m].real*(q[m].x.real*w + x*dwdr*(q[m].value.real - Y[m].real));
					p->fp.x -= c*Y[m].imag*(q[m].x.imag*w + x*dwdr*(q[m].value.imag - Y[m].imag));
					p->fp.y -= c*Y[m].real*(q[m].y.real*w + y*dwdr*(q[m].value.real - Y[m].real));
					p->fp.y -= c*Y[m].imag*(q[m].y.imag*w + y*dwdr*(q[m].value.imag - Y[m].imag));
					p->fp.z -= c*Y[m].real*(q[m].z.real*w + z*dwdr*(q[m].value.real - Y[m].real));
					p->fp.z -= c*Y[m].imag*(q[m].z.imag*w + z*dwdr*(q[m].value.imag - Y[m].imag));
				}
			}
			pij = pij->ilink;
		}
	}
}

void orderSH(SYSTEM*sys, ORDERSH_PARMS*parms, ETYPE*e)
{
	NBR *nbr;
	int i;
	nbr = sys->neighbor;
	particles = nbr->particles;
	L = parms->sph[0]->L;
/*
	vol = sys->box->volume/sys->nglobal;
	scale = cbrt(vol/parms->Vo);
	parms->r1 = parms->r1o*scale;
	parms->r2 = parms->r2o*scale;
*/
	IndexR2 = IndexR1 = -1;
	for (i = 0; i < nbr->nc; i++)
	{
		if (fabs(nbr->rcut[i].value - parms->r2) < 1e-8) IndexR2 = i;
		if (fabs(nbr->rcut[i].value - parms->r1) < 1e-8) IndexR1 = i;
	}
	if (IndexR1 == -1 || IndexR2 == -1)
	{
		printf("Error in order. Index order\n");
		exit(1);
	}
	if (TEST0(sys->loop, parms->globalevalrate))
	{
		e->eion += orderPass1(sys, parms);
		orderPass2(parms, sys->nlocal, sys->nglobal);
	}
	if (TEST0(sys->loop, parms->localevalrate))
	{
		orderSHlocal(sys, parms, e);
	}
}

void orderAccum(COMPLEX*a, COMPLEX*b, double w, int L)
{
	int m;
	for (m = 0; m <= L; m++)
	{
		a[m].real += w*b[m].real;
		a[m].imag += w*b[m].imag;
	}
}

void orderAccumY(COMPLEX*a, YLM*b, double w, int L)
{
	int m;
	for (m = 0; m <= L; m++)
	{
		a[m].real += w*b[m].value.real;
		a[m].imag += w*b[m].value.imag;
	}
}

void orderScale(COMPLEX*a, double s, int L)
{
	int m;
	for (m = 0; m <= L; m++)
	{
		a[m].real *= s;
		a[m].imag *= s;
	}
}

double orderDot(COMPLEX*a, COMPLEX*b, int L)
{
	int m;
	double dot;
	dot = 0.0;
	for (m = 1; m <= L; m++)
	{
		dot += a[m].real*b[m].real + a[m].imag*b[m].imag;
	}
	dot *= 2.0;
	m = 0;
	dot += a[m].real*b[m].real;
	return dot;
}

/**
 *  WARNING: Obtaining memory for static variables such as Q, C, qAccum,
 *  etc from the scratch heap is unsafe since those variable need to be
 *  preserved so that outside actors such as collection_writeBXYZ can
 *  access them.  By convention (as of April 2008) you can't count on
 *  any data staying on the heap once you leave scope.
*/
void orderSHlocal(SYSTEM*sys, ORDERSH_PARMS*parms, ETYPE*e)
{
//	double dot, x, y, z, w, dwdr;
/* 	unsigned qlocal_blk,qnorm_blk,qAccum_blk,Q_blk,Wlocal_blk,C_blk;  */
//	RPAIR *p;
//	PAIRS *pij; 
//	int i, j, k;

	for (int k = 0; k < parms->nL; k++)
	{
		L = parms->sph[k]->L;
		qlocal[k] = (COMPLEX *) ExpandBuffers((void *)qlocal[k], sizeof(COMPLEX)*(L + 1), sys->nion, 1024, LOCATION("orderSHlocal"), "qlocal[0]");
		qnorm[k] = (double *)ExpandBuffers((void *)qnorm[k], sizeof(double), sys->nion, 1024, LOCATION("orderSHlocal"), "qnorm[0]");
	}
	L = parms->sph[0]->L;

	qAccum = (COMPLEX *) ExpandBuffers((void *)qAccum, sizeof(COMPLEX)*(L + 1), sys->nion, 1024, LOCATION("orderSHlocal"), "qAccum");
	Q = (double *)ExpandBuffers((void *)Q, sizeof(double), sys->nion, 1024, LOCATION("orderSHlocal"), "Q");
	Wlocal = (double *)ExpandBuffers((void *)Wlocal, sizeof(double), sys->nion, 1024, LOCATION("orderSHlocal"), "Wlocal");
	C = (int *)ExpandBuffers((void *)C, sizeof(int), sys->nion, 1024, LOCATION("orderSHlocal"), "C");

	for (int k = 0; k < parms->nL; k++)
	{
		L = parms->sph[k]->L;
		for (unsigned i = 0; i < sys->nion*(L + 1); i++)
		{
			qlocal[k][i].real = 0.0;
			qlocal[k][i].imag = 0.0;
		}
	}
	for (unsigned i = 0; i < sys->nion; i++) Q[i] = Wlocal[i] = 0.0;
	for (unsigned i = 0; i < sys->nion; i++) C[i] = 0;

	for (unsigned i = 0; i < sys->nion; i++)
	{
		PAIRS* pij = particles[i].ifirst[IndexR2];
		while (pij  != NULL)
		{
			RPAIR* p = pij->p;
			{
	/*			i = pij->i; */
				int j = pij->j;
				double dummy;
				double w = wfunc(parms, p->r, &dummy);
				double x = p->x;
				double y = p->y;
				double z = p->z;
				for (int k = 0; k < parms->nL; k++)
				{
					L = parms->sph[k]->L;
					sph(parms->sph[k], x, y, z, q);
					orderAccumY(qlocal[k] + (L + 1)*i, q, w, L);
					orderAccumY(qlocal[k] + (L + 1)*j, q, w, L);
				}
			}
			pij = pij->ilink;
		}
	}
	for (unsigned i = 0; i < sys->nion; i++)
	{
		for (int k = 0; k < parms->nL; k++)
		{
			L = parms->sph[k]->L;
			double dot = orderDot(qlocal[k] + (L + 1)*i, qlocal[k] + (L + 1)*i, L);
			dot = sqrt(dot);
			orderScale(qlocal[k] + (L + 1)*i, 1.0/dot, L);
			qnorm[k][i] = dot*sqrt(4.0*M_PI/(2.0*L + 1.0));
		}
	}
	{// scope of k.  I don't know why we didn't just say zero everywhere.
	int k = 0;
	L = parms->sph[k]->L;
/*	for (i = 0; i < sys->nlocal; i++) */  /*Grid BUG*/
	for (unsigned i = 0; i < sys->nion; i++)
	{
		PAIRS* pij = particles[i].ifirst[IndexR2];
		while (pij  != NULL)
		{
			RPAIR* p = pij->p;
			{
				double dummy;
				double w = wfunc(parms, p->r, &dummy);
/*				i=pij->i; */
				int j = pij->j;
				double dot = orderDot(qlocal[k] + (L + 1)*i, qlocal[k] + (L + 1)*j, L);
				Q[i] += dot*w;
				Q[j] += dot*w;
				Wlocal[i] += w;
				Wlocal[j] += w;
				if (dot*w > 0.5)
				{
					C[i] += 1;
					C[j] += 1;
				}
			}
			pij = pij->ilink;
		}
	}
	} // end scope of k
	for (unsigned i = 0; i < sys->nlocal; i++)
	{
		Q[i] /= Wlocal[i];
		for (int k = 0; k < parms->nL; k++) qnorm[k][i] /= Wlocal[i];
	}
	G = (int *)ExpandBuffers((void *)G, sizeof(int), sys->nlocal, 1024, LOCATION("orderSHlocal"), "G");
	for (unsigned i = 0; i < sys->nlocal; i++) G[i] = NOGROUP;
/*
	orderBin(sys, parms);
	orderWrite(sys, parms);
*/
/*	orderCluster(sys, parms);*/
}

#ifdef COMPILE_UNUSED 
static double Qc[] = { -0.5, 0.75, 0.87, 0.95, 1.0 };
void orderWrite(SYSTEM*sys, ORDERSH_PARMS*parms)
{
	int id = getRank(0);
	char filename[1024];
	sprintf(filename, "Axyz.%6.6d", id);
	FILE* Axyzfile = fopen(filename, "w");
	int nL = parms->nL;
	double* rx = sys->collection->state->rx;
	double* ry = sys->collection->state->ry;
	double* rz = sys->collection->state->rz;
	gid_type* label = sys->collection->state->label;
	for (unsigned j = 0; j < sys->nlocal; j++)
	{
		if (Q[j] > Qc[CRYSTAL])
		{
			float xs = rx[j];
			float ys = ry[j];
			float zs = rz[j];
			float Qs = Q[j];
			fprintf(Axyzfile, "%12llu %10.4f %10.4f %10.4f %10.4f %2d %2d", label[j], xs, ys, zs, Qs, C[j], nL);
			for (int k = 0; k < nL; k++)
			{
				float qs = qnorm[k][j];
				fprintf(Axyzfile, " %8.3f", qs);
			}
			fprintf(Axyzfile, "\n");
		}
	}
}
void orderBin(SYSTEM*sys, ORDERSH_PARMS*parms)
{
	int id = getRank(0);
	nL = parms->nL;
	int Ncount[3];
	Ncount[0] = Ncount[1] = Ncount[2] = 0;
	const double deltaBin = 0.02;
	const double deltaBinI = 1.0/deltaBin;
	const int nbin = 80;
	double* bindata = (double *)ExpandBuffers((void *)bindata, sizeof(double), nbin*(2 + 3*nL), 1024, LOCATION("orderBin"), "bindata");
	double* binR = NULL;
	if (id == 0) binR = (double *)ExpandBuffers((void *)binR, sizeof(double), nbin*(2 + 3*nL), 1024, LOCATION("orderBin"), "binR");
	for (int k = 0; k < nbin; k++) for (int j = 0; j < nL + 2; j++) bindata[k + j*nbin] = 0.0;
	for (unsigned i = 0; i < sys->nlocal; i++)
	{
		unsigned kType = LIQUID;
		if (Q[i] > Qc[INTERFACE]) kType = INTERFACE;
		if (Q[i] > Qc[CRYSTAL]) kType = CRYSTAL;
		Ncount[kType]++;
		for (int j = 0; j < nL; j++)
		{
			int ibin = deltaBinI*(qnorm[j][i]);
			if (ibin > 79) ibin = 79;
			ibin += nbin*(j + nL*kType);
			bindata[ibin] += 1.0;
		}
		int ibin = 5*C[i];
		if (ibin > 79) ibin = 79;
		ibin += nbin*(3*nL);
		bindata[ibin] += 1.0;

		ibin = deltaBinI*(Q[i] + 0.5);
		if (ibin < 0) ibin = 0;
		ibin += nbin*(3*nL + 1);
		bindata[ibin] += 1.0;
		if (kType < CRYSTAL) G[i] = kType;
	}
	MPI_Reduce(bindata, binR, nbin*(2 + 3*nL), MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
	printf("Norms %d %d %d %d\n", Ncount[LIQUID], Ncount[INTERFACE], Ncount[CRYSTAL], sys->nlocal);
	fflush(stdout);
/*
	{
		sprintf(filename,"bin#%6.6d",getRank(0));
		file = fopen(filename,"w");
		for (i=0;i<nbin;i++)
		{
			fprintf(file,"%4.2f",(i+0.5)*deltaBin);
			for (k=0;k<3;k++) for (j=0;j<nL+1;j++) fprintf(file," %6.0f",bindata[i+nbin*(j+k*nL)]);
			fprintf(file," %5.2f",(i+0.5)*deltaBin-0.5);
			fprintf(file," %6.0f",bindata[i+nbin*(3*nL+1)]);
			fprintf(file,"\n");
		}
		fclose(file);
	}
*/
	if (id == 0)
	{
		FILE* file = fopen("distqQ", "w");
		for (int i = 0; i < nbin; i++)
		{
			fprintf(file, "%4.2f", (i + 0.5)*deltaBin);
			for (int k = 0; k < 3; k++) for (int j = 0; j < nL + 1; j++) fprintf(file, " %6.0f", binR[i + nbin*(j + k*nL)]);
			fprintf(file, " %5.2f", (i + 0.5)*deltaBin - 0.5);
			fprintf(file, " %6.0f", binR[i + nbin*(3*nL + 1)]);
			fprintf(file, "\n");
		}
		fclose(file);
	}
}

void orderCluster(SYSTEM*sys, ORDERSH_PARMS*parms)
{
	COMPLEX* qAve = Y;
	int id = getRank(0);
/*
	sprintf(filename,"bxyz.%6.6d",id);
	bxyzfile=fopen(filename,"w");
*/
	char filename[1024];
	sprintf(filename, "cluster.%6.6d", id);
	FILE* clusterfile = fopen(filename, "w");
	nL = parms->nL;
	const int k0 = 0;
	L = parms->sph[k0]->L;
	for (unsigned i = 0; i < sys->nlocal; i++)
		for (int j = 0; j < L + 1; j++)
			qAccum[(L + 1)*i + j] = qlocal[k0][(L + 1)*i + j];
	for (unsigned i = 0; i < sys->nlocal; i++)
	{
		PAIRS* pij = particles[i].ifirst[IndexR1];
		while (pij  != NULL)
		{
		/*	i = pij->i; */
			int j = pij->j;
			double dot = orderDot(qlocal[k0] + (L + 1)*i, qlocal[k0] + (L + 1)*j, L);
			double w = 1.0;
			if (dot > 0.95 && Q[i] > Qc[HIGHORDER] && Q[j] > Qc[HIGHORDER])
			{
				orderAccum(qAccum + (L + 1)*i, qlocal[k0] + (L + 1)*j, w, L);
				orderAccum(qAccum + (L + 1)*j, qlocal[k0] + (L + 1)*i, w, L);
			}
			pij = pij->ilink;
		}
	}
	for (unsigned i = 0; i < sys->nlocal; i++)
	{
		double dot = sqrt(orderDot(qAccum + (L + 1)*i, qAccum + (L + 1)*i, L));
		orderScale(qAccum + (L + 1)*i, 1.0/dot, L);
	}
	int ngroup = 0;
	int NT = 0;
	double* rx = sys->collection->state->rx;
	double* ry = sys->collection->state->ry;
	double* rz = sys->collection->state->rz;
	gid_type* label = sys->collection->state->label;
	unsigned i;
	for (i = 0; i < sys->nlocal; i++) if (Q[i] > Qc[HIGHORDER]) break;
	while (ngroup < 64 && i < sys->nlocal)
	{
		int n = 0;
//unused		int group = getRank(0)*1000 + ngroup;
		for (int m = 0; m < L + 1; m++) qAve[m].real = qAve[m].imag = 0.0;
		for (unsigned j = 0; j < sys->nlocal; j++)
		{
			if (Q[j] > Qc[HIGHORDER] && G[j] == NOGROUP)
			{
				double dot = orderDot(qAccum + (L + 1)*i, qlocal[k0] + (L + 1)*j, L);
				if (dot > 0.95) orderAccum(qAve, qlocal[k0] + (L + 1)*j, 1.0, L);
			}
		}
		double dot = sqrt(orderDot(qAve, qAve, L));
		orderScale(qAve, 1.0/dot, L);
		THREE_VECTOR r;
		double r2;
		r.x = r.y = r.z = r2 = 0.0;
// unused		int ptr = ftell(clusterfile);
		for (unsigned j = 0; j < sys->nlocal; j++)
		{
			if (Q[j] > Qc[CRYSTAL] && G[j] == NOGROUP)
			{
				dot = orderDot(qAve, qlocal[k0] + (L + 1)*j, L);
				if (dot > 0.85)
				{
					G[j] = ngroup + 100;
					double x = rx[j] - rx[i];
					double y = ry[j] - ry[i];
					double z = rz[j] - rz[i];
					nearestImage(&x, &y, &z);
					r.x += x;
					r.y += y;
					r.z += z;
					r2 += x*x + y*y + z*z;
					float xs = rx[j];
					float ys = ry[j];
					float zs = rz[j];
					float Qs = Q[j];
					float dots = dot;
					fprintf(clusterfile, "%llu %d %f %f %f %f %f %d %d", label[j], ngroup, xs, ys, zs, dots, Qs, C[j], nL);
					for (int k = 0; k < nL; k++)
					{
						float qs = qnorm[k][j];
						fprintf(clusterfile, " %f", qs);
					}
					fprintf(clusterfile, "\n");
					n++;
				}
			}
		}
		r.x /= n;
		r.y /= n;
		r.z /= n;
		r2 /= n;
		double rmsR = sqrt(n*(r2 - r.x*r.x - r.y*r.y - r.z*r.z)/(n - 1));
		r.x += rx[i];
		r.y += ry[i];
		r.z += rz[i];
		cluster[ngroup].label = ngroup;
		cluster[ngroup].size = n;
		cluster[ngroup].radius = 0.0;
		cluster[ngroup].Rrms = rmsR;
		cluster[ngroup].Rave = r;
		for (int m = 0; m < L + 1; m++)
			cluster[ngroup].qAve[m] = qAve[m];
		ngroup++;
		NT += n;
		for (; i < sys->nlocal; i++)
			if (Q[i] > Qc[HIGHORDER] && G[i] == NOGROUP) break;
	}
	for (unsigned j = 0; j < sys->nlocal; j++)
	{
		int k = G[j] - 100;
		if (k >= 0)
		{
			double x = rx[j] - cluster[k].Rave.x;
			double y = ry[j] - cluster[k].Rave.y;
			double z = rz[j] - cluster[k].Rave.z;
			nearestImage(&x, &y, &z);
			double r2 = x*x + y*y + z*z;
			if (r2 > cluster[k].radius) cluster[k].radius = r2;
		}
	}
	for (int k = 0; k < ngroup; k++) cluster[k].radius = sqrt(cluster[k].radius);
	fclose(clusterfile);
	{
		static int *displs = NULL;
		static int* ngroups = NULL;
		if (id == 0)
		{
			ngroups = (int *)ExpandBuffers((void *)ngroups, sizeof(int), getSize(0), 1024, LOCATION("orderCluster"), "ngroups");
			displs = (int *)ExpandBuffers((void *)displs, sizeof(int), getSize(0), 1024, LOCATION("orderCluster"), "displs");
		}
		MPI_Gather(&ngroup, 1, MPI_INT, ngroups, 1, MPI_INT, 0, COMM_LOCAL);
		if (id == 0)
		{
			displs[0] = 0;
			for (int j = 1; j < getSize(0); j++)
				displs[j] = displs[j - 1] + ngroups[j - 1];
			ngroup = 0;
		}
		MPI_Gatherv(cluster, ngroup, CLUSTER_TYPE, cluster, ngroups, displs, CLUSTER_TYPE, 0, COMM_LOCAL);
		if (id == 0)
		{
			clusterfile = fopen("clusterinfo", "w");
			for (int j = 0; j < getSize(0); j++)
				ngroup += ngroups[j];
			printf("%d ", ngroup);
			ngroup = clusterMerge(ngroup, L, cluster);
			printf("%d ", ngroup);
			for (int j = 0; j < ngroup; j++)
			{
				double r = cluster[j].Rrms;
				double vol = 4*M_PI/3.0*r*r*r/cluster[j].size;
				fprintf(clusterfile, "%4d %4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f",
					cluster[j].label, cluster[j].size, cluster[j].radius, cluster[j].Rrms, vol, cluster[j].Rave.x, cluster[j].Rave.y, cluster[j].Rave.z);
				for (int m = 0; m < L + 1; m++)
					fprintf(clusterfile, " %7.4f %7.4f", cluster[j].qAve[m].real, cluster[j].qAve[m].imag);
				fprintf(clusterfile, "\n");
				/*fprintf(groupfile," %s %d %d %d\n",filename,1,n,ptr); */
			}
			fclose(clusterfile);
		}
	}
	exit(0);
}

int clusterMerge(int n, int L, CLUSTER*cluster)
{
	double dot, doti, dotj;
	COMPLEX *si, *sj;
	double x, y, z, r2;
	int i, j, k, ni, nj, m;
	for (i = 0; i < n; i++)
	{
		si = cluster[i].qAve;
		for (j = i + 1; j < n; j++)
		{
			sj = cluster[j].qAve;
			doti = orderDot(si, si, L);
			dotj = orderDot(sj, sj, L);
			dot = orderDot(si, sj, L)/sqrt(doti*dotj);
			x = cluster[j].Rave.x - cluster[i].Rave.x;
			y = cluster[j].Rave.y - cluster[i].Rave.y;
			z = cluster[j].Rave.z - cluster[i].Rave.z;
			nearestImage(&x, &y, &z);
			r2 = x*x + y*y + z*z;
			if (dot > 0.85)
			{
				ni = cluster[i].size;
				nj = cluster[j].size;
				printf("dot = %f r= %f Merging %d %d %d %d\n", dot, sqrt(r2), i, j, ni, nj);
/*
				printf(" %4d %4d %10.4f %10.4f %10.4f %10.4f %10.4f\n", 				cluster[i].label,cluster[i].size,cluster[i].radius,cluster[i].Rrms,cluster[i].Rave.x,cluster[i].Rave.y,cluster[i].Rave.z);
				printf(" %4d %4d %10.4f %10.4f %10.4f %10.4f %10.4f\n", 				cluster[j].label,cluster[j].size,cluster[j].radius,cluster[j].Rrms,cluster[j].Rave.x,cluster[j].Rave.y,cluster[j].Rave.z);
*/
				for (m = 0; m < L + 1; m++)
				{
					cluster[i].qAve[m].real = (ni*si[m].real + nj*sj[m].imag)/(ni + nj);
					cluster[i].qAve[m].imag = (ni*si[m].imag + nj*sj[m].imag)/(ni + nj);
				}
/*
				vi = cluster[i].Rrms;
				vj = cluster[j].Rrms;
				cluster[i].Rrms = sqrt((vi*vi*ni + vj*vj*nj + ni*nj*r2/(ni+nj))/(ni+nj));
				cluster[i].Rave.x += nj*x /(ni+nj);
				cluster[i].Rave.y += nj*y /(ni+nj);
				cluster[i].Rave.z += nj*z /(ni+nj);
*/
				cluster[i].size = (ni + nj);
				for (k = j; k < n - 1; k++)
					cluster[k] = cluster[k + 1];
				n--;
			}
		}
	}
	return n;
}
#endif

int orderGetnL(void)
{
	return nL;
}

double *orderGetQ(void)
{
	return Q;
}
int *orderGetLv(void)
{
	return Lv;
}

double **orderGetqnorm(void)
{
	return qnorm;
}

COMPLEX *orderGetqAccum(void)
{
	return NULL;
	return qAccum;
}

int *orderGetC(void)
{
	return C;
}
void writeqlocal(SIMULATE*simulate)
{
/* 	PFILE *file; */
/* 	char filename[512],field_names[1024],field_types[1024],name[16]; */
/* 	int i, l,L,k, m, n, nfields;  */
/* 	unsigned  lrec,nlocal; */
/* 	gid_type  nglobal; */
/* 	float b[1024]; */
	gid_type nglobal = simulate->system->nglobal;
	unsigned nlocal = simulate->system->nlocal;
	for (int k=0;k<nL;k++)
	{
		int L=Lv[k];
		char filename[512],field_names[1024],field_types[1024],name[16];
		snprintf(filename, 512,"%s/q%d", simulate->snapshotdir,L);
		snprintf(name, 16,"q%d", L);
		PFILE* file = Popen(filename, "w", COMM_LOCAL);
		snprintf(field_names,1024,"checksum ");
		snprintf(field_types,1024,"u4 ");
		for (int m=0;m<=L;m++)
		{
			int l=strlen(field_names);
			snprintf(field_names+l,1023-l,"q%dr[%d] q%di[%d] ",L,m,L,m);
			l=strlen(field_types);
			snprintf(field_types+l,1023-l,"f4 f4 ");
		}
		int nfields = 1 + 2*(L+1); 
/*
		lrec =1; 
		while (lrec<4*nfields) lrec *= 2; 
*/
		unsigned lrec = 4*nfields; 
		PioSet(file, "recordLength", lrec);
		PioSet(file, "datatype", FIXRECORDBINARY);
		PioSet(file, "numberRecords", nglobal);
		PioSet(file, "checksum", CRC32);
		PioSet(file, "nfields", nfields);
		PioSet(file, "field_names", field_names);
		PioSet(file, "field_types", field_types);
		PioReserve(file, lrec*nlocal + 2048);
		if (getRank(0) == 0) write_fileheader(file, simulate, name );
		for (unsigned i = 0; i < nlocal; i++)
		{
			float b[1024];
			unsigned n=1;
			for (int m=0;m<=L;m++)
			{
				b[n++] = qlocal[k][(L+1)*i+m].real;
				b[n++] = qlocal[k][(L+1)*i+m].imag;
			}
			for (;4*n<lrec;n++) b[n]=0;
   			unsigned checksum = checksum_crc32_table((unsigned char *)b+1,lrec-4);
      		copyBytes(b, &checksum, 4);
			Pwrite(b, lrec, 1, (PFILE *)file);
		}
		Pclose(file);
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
