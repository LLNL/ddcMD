#include "md2ddc.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "three_algebra.h"
#include "ddc.h"
#include "species.h"
#include "group.h"
static int *nlocalptr = NULL;
static int *nionptr = NULL;
static double *rxptr = NULL, *ryptr = NULL, *rzptr = NULL;
static double *vxptr = NULL, *vyptr = NULL, *vzptr = NULL;
static double *fxptr = NULL, *fyptr = NULL, *fzptr = NULL;
static double *energyptr = NULL, *virialptr = NULL;
static THREE_SMATRIX *sionptr = NULL; 
static int *typeptr = NULL;
static gid_type  *labelptr = NULL;
static SPECIES **speciesptr = NULL;
static GROUP **groupptr = NULL;
SPECIES *species_by_index(SPECIES **, int);
GROUP *group_by_index(GROUP **, int);
void setparticleptr_(double *rx, double *ry, double *rz, int *type, SPECIES ** species, gid_type *label, int *n)
{
	rxptr = rx;
	ryptr = ry;
	rzptr = rz;
	typeptr = type;
	labelptr = label;
	speciesptr = species;
	nlocalptr = n;
}
int get_ParticleNumber()
{
	return *nlocalptr;
}

void get_particles(PARTICLE*ddc_particle, int ddc_n)
{
	int i;
	for (i = 0; i < ddc_n; i++)
	{
		ddc_particle[i].r.x = rxptr[i];
		ddc_particle[i].r.y = ryptr[i];
		ddc_particle[i].r.z = rzptr[i];
		ddc_particle[i].type = typeptr[i];
		ddc_particle[i].global_index = labelptr[i];
	}
}

void put_particles(PARTICLE*ddc_particle, int nlocal, int nremote)
{
	for (int i = 0; i < nlocal + nremote; i++)
	{
		rxptr[i] = ddc_particle[i].r.x;
		ryptr[i] = ddc_particle[i].r.y;
		rzptr[i] = ddc_particle[i].r.z;
		int index = typeptr[i] = ddc_particle[i].type;
		speciesptr[i] = species_by_index(NULL, index & 0xffff);
		groupptr[i] = group_by_index(NULL, index >> 16);
		labelptr[i] = ddc_particle[i].global_index;
	}
}

void *GetRecvLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE)
{
	int blockcounts[3];
	MPI_Datatype type[3];
	MPI_Aint displs[3];
	type[0] = MPI_DOUBLE;
	type[1] = MPI_DOUBLE;
	type[2] = MPI_DOUBLE;
	blockcounts[0] = nRecv;
	blockcounts[1] = nRecv;
	blockcounts[2] = nRecv;
	MPI_Get_address(rxptr + RecvStart, &displs[0]);
	MPI_Get_address(ryptr + RecvStart, &displs[1]);
	MPI_Get_address(rzptr + RecvStart, &displs[2]);
	int n =3; 
	for (int i = n-1; i >= 0; i--) displs[i] -= displs[0];
	if (*(RECVTYPE) != 0) MPI_Type_free(RECVTYPE);
	MPI_Type_create_struct(n, blockcounts, displs, type, RECVTYPE);
	MPI_Type_commit(RECVTYPE);
	return rxptr + RecvStart;
}

void *GetForceLocation(int n, int Start, MPI_Datatype*TYPE,unsigned pCalculate)
{
   
	MPI_Datatype type[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	int  blockcounts[] ={n,n,n,n,n,0};
	MPI_Aint displs[6];
   int m=0;
	MPI_Get_address(fxptr + Start, &displs[m++]);
	MPI_Get_address(fyptr + Start, &displs[m++]);
	MPI_Get_address(fzptr + Start, &displs[m++]);
	if (pCalculate & 1) MPI_Get_address(energyptr + Start, &displs[m++]);
	if (pCalculate & 2) MPI_Get_address(virialptr + Start, &displs[m++]);
	if (pCalculate & 4) 
   {
      blockcounts[m]=6*n;
      MPI_Get_address(sionptr + Start, &displs[m++]);
   }
	for (int i = m-1; i >= 0; i--) displs[i] -= displs[0];
	if (*(TYPE) != 0) MPI_Type_free(TYPE);
	MPI_Type_create_struct(m, blockcounts, displs, type, TYPE);
	MPI_Type_commit(TYPE);
	return fxptr + Start;
}

void zeroremote(int nlocal, int nremote)
{
	int i;
	for (i = 0; i < nremote; i++)
	{
		rxptr[i + nlocal] = 0;
		ryptr[i + nlocal] = 0;
		rzptr[i + nlocal] = 0;
	}
}
void fillsendbuf(double *send, int nlist, int *list)
{
	double *x, *y, *z;
	int i;
	x = send;
	y = x + nlist;
	z = y + nlist;
	for (i = 0; i < nlist; i++)
	{
		x[i] = rxptr[list[i]];
		y[i] = ryptr[list[i]];
		z[i] = rzptr[list[i]];
	}
}
void *GetVelocityLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE)
{
	int i, blockcounts[3];
	MPI_Datatype type[3];
	MPI_Aint displs[3];
	type[0] = MPI_DOUBLE;
	type[1] = MPI_DOUBLE;
	type[2] = MPI_DOUBLE;
	blockcounts[0] = nRecv;
	blockcounts[1] = nRecv;
	blockcounts[2] = nRecv;
	MPI_Get_address(vxptr + RecvStart, &displs[0]);
	MPI_Get_address(vyptr + RecvStart, &displs[1]);
	MPI_Get_address(vzptr + RecvStart, &displs[2]);
	int n =3; 
	for (i = n-1; i >= 0; i--) displs[i] -= displs[0];
	if (*(RECVTYPE) != 0) MPI_Type_free(RECVTYPE);
	MPI_Type_create_struct(n, blockcounts, displs, type, RECVTYPE);
	MPI_Type_commit(RECVTYPE);
	return vxptr + RecvStart;
}
void Velocityfillsendbuf(double *send, int nlist, int *list)
{
	double *vx,*vy,*vz;
	int i;
	vx = send;
	vy = vx + nlist;
	vz = vy + nlist;
	for (i = 0; i < nlist; i++)
	{
		vx[i] = vxptr[list[i]];
		vy[i] = vyptr[list[i]];
		vz[i] = vzptr[list[i]];
	}
}
void forceaccum0(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
	}
}
void forceaccum1(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *e = z + nlist;
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		energyptr[list[i]] += e[i];
	}
}
void forceaccum2(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *v = z + nlist;
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		virialptr[list[i]] += v[i];
	}
}
void forceaccum3(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *e = z + nlist;
	double *v = e + nlist;
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		energyptr[list[i]] += e[i];
		virialptr[list[i]] += v[i];
	}
}
void forceaccum4(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	THREE_SMATRIX *s = (THREE_SMATRIX *)(z + nlist);
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		sionptr[list[i]].xx += s[i].xx;
		sionptr[list[i]].xy += s[i].xy;
		sionptr[list[i]].xz += s[i].xz;
		sionptr[list[i]].yy += s[i].yy;
		sionptr[list[i]].yz += s[i].yz;
		sionptr[list[i]].zz += s[i].zz;
	}
}
void forceaccum5(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *e = z + nlist;
	THREE_SMATRIX *s = (THREE_SMATRIX *)(e + nlist);
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		energyptr[list[i]] += e[i];
		sionptr[list[i]].xx += s[i].xx;
		sionptr[list[i]].xy += s[i].xy;
		sionptr[list[i]].xz += s[i].xz;
		sionptr[list[i]].yy += s[i].yy;
		sionptr[list[i]].yz += s[i].yz;
		sionptr[list[i]].zz += s[i].zz;
	}
}
void forceaccum6(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *v = z + nlist;
	THREE_SMATRIX *s = (THREE_SMATRIX *)(v + nlist);
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		virialptr[list[i]] += v[i];
		sionptr[list[i]].xx += s[i].xx;
		sionptr[list[i]].xy += s[i].xy;
		sionptr[list[i]].xz += s[i].xz;
		sionptr[list[i]].yy += s[i].yy;
		sionptr[list[i]].yz += s[i].yz;
		sionptr[list[i]].zz += s[i].zz;
	}
}
void forceaccum7(double *recv, int nlist, int *list)
{
	double *x = recv;
	double *y = x + nlist;
	double *z = y + nlist;
	double *e = z + nlist;
	double *v = e + nlist;
	THREE_SMATRIX *s = (THREE_SMATRIX *)(v + nlist);
	for (int i = 0; i < nlist; i++)
	{
		fxptr[list[i]] += x[i];
		fyptr[list[i]] += y[i];
		fzptr[list[i]] += z[i];
		energyptr[list[i]] += e[i];
		virialptr[list[i]] += v[i];
		sionptr[list[i]].xx += s[i].xx;
		sionptr[list[i]].xy += s[i].xy;
		sionptr[list[i]].xz += s[i].xz;
		sionptr[list[i]].yy += s[i].yy;
		sionptr[list[i]].yz += s[i].yz;
		sionptr[list[i]].zz += s[i].zz;
	}
}

void writemg(FILE*file, int start, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		fprintf(file, "writemg: %u %"PRIu64" %f %f %f\n", i, labelptr[i + start], rxptr[i + start], ryptr[i + start], rzptr[i + start]);
	}
	fflush(file);
}

void mgpt2ddcmemset(STATE*state)
{
	rxptr = state->rx;
	ryptr = state->ry;
	rzptr = state->rz;
	vxptr = state->vx;
	vyptr = state->vy;
	vzptr = state->vz;
	labelptr = state->label;
	typeptr = state->atomtype;
	speciesptr = state->species;
	groupptr = state->group;
	fxptr = state->fx;
	fyptr = state->fy;
	fzptr = state->fz;
	energyptr = state->potentialEnergy;
	virialptr = state->virial;
	sionptr = state->sion;
	nlocalptr = &state->nlocal; 
	nionptr = &state->nion; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
