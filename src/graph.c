#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "three_algebra.h"
#include "mgpt.h"
#include "simulate.h"
#include "printinfo.h"
#include "system.h"
#include "units.h"
#include "mpiUtils.h"
#include "format.h"

enum GRAPHS {nlocal_, nremote_, nion_, npair_, nRpair_, nps_, npsllMAX_, npsMAX_, n4_, n43_, n33_, n2_}; 
static double flopcnt; 
unsigned *mgptgetNgraphs(void);

#ifdef COMPILE_UNUSED
static double getflopcnt(void) 
{
	return  flopcnt; 
}
#endif
void graphWrite(SIMULATE*simulate, FILE*file)
{
	SYSTEM *sys;
   if (simulate->printinfo->printGraphs == 0) return; 
	int ibuf[64], minbuf[64], maxbuf[64];
	static int flops[12],cdiv=4,csqrt=12,order;
	LONG64 sendbuf[64], sumbuf[64];
	double avebuf[64];
	static int header = 1;
	int i, j;
	double time;
	unsigned id, size, nlocal, nremote, nion, npair, nRpair, *ngraphs;
	id = getRank(0);
	size = getSize(0);
	SIGNED64 loop = simulate->system->loop;
	time = units_convert(simulate->time,NULL,"t");
	LONG64 nglobal = simulate->system->nglobal;
	nlocal = simulate->system->nlocal;
	nremote = simulate->system->nremote;
	nion = simulate->system->nion;
	sys = simulate->system;
	npair = simulate->system->neighbor->npairs;
	nRpair = simulate->system->neighbor->nRpair;
/*
		np0=n2;
			c4 =  (2*250 + 3 * 406 + 450  + 5)*n4 ;   
			c33 = (2*250 + 406+3)*n33 ;
			chamltn = 522*nps;
			c4 =  (2*686 + 3 * 986 + 986 + 5)*n4 ; 
			c33 =  (2*986 + 989+3)*n33 ;
			chamltn = 2623*nps;
			c43 =  73*n43;
			c2  =  (19 + cdiv)*n2;
			cpairlist =   (11 + csqrt)*8*npair*f;
			cpairUpdate = (8+csqrt)*npair*(1-f);
			cfker=21 *np0 ; 
			cfsumX=12*nRpair;
			cnt = c4+c43+c33+c2+cpairlist+cpairUpdate+cfker+cfsumX+chamltn;
			cnt /= nlocal;
*/
	ngraphs = mgptgetNgraphs();
	ibuf[0] = sendbuf[0] = nlocal;
	ibuf[1] = sendbuf[1] = nremote;
	ibuf[2] = sendbuf[2] = nion;
	ibuf[3] = sendbuf[3] = npair;
	ibuf[4] = sendbuf[4] = nRpair;
	for (i = 0; i < 7; i++) ibuf[i + 5] = sendbuf[i + 5] = ngraphs[i];
	MPI_Reduce(ibuf, minbuf, 12, MPI_INT, MPI_MIN, 0, COMM_LOCAL);
	MPI_Reduce(ibuf, maxbuf, 12, MPI_INT, MPI_MAX, 0, COMM_LOCAL);
	MPI_Reduce(sendbuf, sumbuf, 12, MPI_LONG_LONG_INT, MPI_SUM, 0, COMM_LOCAL);
/*
	static LONG64 HWLast=0;
	LONG64 HW = getHPM_FlopCnt("LOOP"); 
	if (id == 0) printf("HW: %e %e %e\n",HW*1.0,HWLast*1.0,(HW-HWLast)*1.0);
	double HWflopcnt = ((double)(HW - HWLast))/simulate->printrate + sumbuf[n4_]*414 + sumbuf[n33_]*228;
	if (id ==0) printf("HWflopcnt=%e %e\n",HW*1.0/(simulate->printrate*nglobal),HWflopcnt/nglobal);
	HWLast = HW; 
*/
	if (id != 0) return;
	for (i = 0; i < 12; i++) avebuf[i] = (1.0*sumbuf[i])/size;
	flops[npsllMAX_] = 0;
	flops[npsMAX_] = 0;
	if (header)
	{
		fprintf(file, "%-12s %12s %12s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s",
			"#       loop", "time", "nglobal", "nlocal", "nremote", "nion", "npair", "nRpair", "nps", "npsllMax", "npsMAX", "n4", "n43", "n33", "n2");
		header = 0;
		for (i = 0; i < sys->npotential; i++)
		{
			if (sys->potential[i] && sys->potential[i]->itype == MGPT)
			{
				double f=0.1; 
				order = ((MGPT_PARMS*)sys->potential[0]->parms)->hamltninfo->order;
				flops[nlocal_] = 0; 
				flops[nremote_] = 0; 
				flops[nion_] = 0; 
				flops[npair_] = (11+csqrt)*8*f+(8+csqrt)*1.0-f; 
				flops[nRpair_] = 18;
				if (order==5)
				{
					flops[nps_] = 522;     
					flops[n4_] = (2*250 + 3 * 406 + 456 + 5);
					flops[n43_] = 73 ; 
					flops[n33_] = (2*406+456+3) ; 
				}
				if (order==7)
				{
					flops[nps_] = 2623;     
					flops[n4_] = (2*686 + 3 * 986 + 1084 + 5) ; 
					flops[n43_] = 73 ; 
					flops[n33_] = (2*986+1084+3) ; 
				}
				flops[n2_] =  22 + cdiv + 22;   /*vfunc  +  fker */
			}
			if (sys->potential[i] && sys->potential[i]->itype == ORDERSH)
			{
			}
		}
		if (simulate->integrator->itype == NPTGLF)
		{
		}
		fprintf(file, "\n");
	}
	flopcnt=0; 
	for (j = 0; j < 12; j++) flopcnt += sumbuf[j]*flops[j];
	fprintf(file, "min ");
	fprintf(file, loopFormat(), loop);
	fprintf(file, " %12.2f %12"PRIu64, time, nglobal);
	for (i = 0; i < 12; i++) fprintf(file, " %8d", minbuf[i]);
	fprintf(file, "\n");
	fprintf(file, "max ");
	fprintf(file, loopFormat(), loop);
	fprintf(file, " %12.2f %12"PRIu64, time, nglobal);
	for (i = 0; i < 12; i++) fprintf(file, " %8d", maxbuf[i]);
	fprintf(file, "\n");
	fprintf(file, "ave ");
	fprintf(file, loopFormat(), loop);
	fprintf(file, " %12.2f %12"PRIu64, time, nglobal);
	for (i = 0; i < 12; i++) fprintf(file, " %8.0f", avebuf[i]);
	fprintf(file, "\n");
/*
	fprintf(file, "fac ");
	fprintf(file, loopFormat(), loop);
	fprintf(file, " %12.2f %12"PRIu64, time, nglobal);
	for (i = 0; i < 12; i++) fprintf(file, " %8.4f", sumbuf[i]*flops[i]*1.0/flopcnt);
*/
/*
	fprintf(file, " %8.0f", flopcnt/nglobal);
	fprintf(file, " %8.0f", HWflopcnt/nglobal);
	printf("flopcnt = %16.8f\n ", flopcnt/nglobal);
	printf("HWflopcnt = %16.8f %f %f\n ", HWflopcnt/nglobal,HWflopcnt/flopcnt,HWflopcnt-flopcnt);
*/
	fprintf(file, "\n\n");
	fflush(file);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
