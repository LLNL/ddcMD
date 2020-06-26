#include "ddc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "geom.h"
#include "object.h"
#include "particle.h"
#include "box.h"
#include "ddcMalloc.h"
#include "neighbor.h"
#include "expandbuffer.h"
#include "mpiUtils.h"
#include "preduce.h"
#define TAG 1


//static char *DOMAINTYPENAMES[]  = {"NONE","CUBIC", "BCC", "FCC", "HEXAGONAL", "NONLATTICE" };




void *GetRecvLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE);

void sortsetparticle(PARTICLE*pvalue);

void ddcFindPairs(DDC*ddc)
{
	sortsetparticle(ddc->particles);
	profile(BOX, START);

	PARTICLESET *particleset = NULL;
	PARTICLE *dp = ddc->particles;
   const THREE_VECTOR* domainCenterP = &(ddc->domains.domains[ddc->domain_id].center);
	int pSize = sizeof(PARTICLE);
	particleset = 
		ParticleSet(particleset, &(dp->r.x), &(dp->r.y), &(dp->r.z), pSize, 
						&(dp->global_index), pSize, &(dp->type), pSize,
						&ddc->number_particles, &(ddc->number_local), domainCenterP);
	*particleset->number_local = ddc->number_local;
	*particleset->number_particles = ddc->number_particles;

	GEOM* geom = geom_init(ddc->rcut, box_getBox(NULL)); 
	
	GeomBox(geom, particleset);
	profile(BOX, END);

	double rcut2 = SQ(ddc->rcut);
	int nRemoteProcessors= ddc->nRemoteProcessors ;
	int nlocal = ddc->number_local;
	int nremote= ddc->number_remote;

	MPI_Request sendRequest[nRemoteProcessors+1];
	int nremoteNew =0; 
	int recvStartNew =nlocal; 
	int keepRecv[nremote+1]; 
	for (int i = 0; i < nremote; i++)keepRecv[i] =0; 
	for (int ii=0; ii<nRemoteProcessors; ++ii)
	{
	    int remoteId=ddc->CommInfo[ii].domain_id ;
		int recvStart = ddc->CommInfo[ii].RecvStart; 
		int nRecv = ddc->CommInfo[ii].nRecv; 

		if (nRecv > 0)
		{
			for (int i = recvStart; i < recvStart+nRecv; i++)
			{
				THREE_VECTOR r = ddc->particles[i].r;
				GEOMBOX *box = geom->pinfo[i].box;
				int j=-1;
				for (int k = 0; k < box->nn; k++)
				{
					GEOMBOX *box_neigh = box + box->nlist[k];
					j = box_neigh->firstlocal;
					while (j != -1)
					{
						double x = r.x - ddc->particles[j].r.x;
						double y = r.y - ddc->particles[j].r.y;
						double z = r.z - ddc->particles[j].r.z;
						double r2 = x*x + y*y + z*z;
						if (r2 > rcut2) nearestImage_fast(&x,&y,&z); r2 = x*x + y*y + z*z;
 						if (r2 < rcut2)  
						{
							keepRecv[i-nlocal] = 1; 
							k = box->nn; 
							break; 
						}
						j = geom->pinfo[j].next;
					}
				}
			}
		}

		MPI_Isend(keepRecv+recvStart-nlocal, nRecv, MPI_INT, remoteId, TAG, COMM_LOCAL, sendRequest + ii);
		int nRecvNew = 0; 
		ddc->CommInfo[ii].RecvStart = recvStartNew ;
		for (int i = 0; i < nRecv; i++) if (keepRecv[i+recvStart-nlocal]) ddc->particles[recvStartNew+nRecvNew++] = ddc->particles[i+recvStart];
		ddc->CommInfo[ii].nRecv = nRecvNew;
		ddc->CommInfo[ii].RecvPtr = GetRecvLocation(nRecvNew, recvStartNew, &(ddc->CommInfo[ii].RecvType));
		nremoteNew += nRecvNew; 
		recvStartNew += nRecvNew; 

	}
	geom_free(geom);
	geom=NULL;
	ddcFree(particleset);
	ddc->nSend = 0;
	int sendStartNew=0; 
   profile(WAIT_FPAIRS, START);
   for (int ii=0; ii<nRemoteProcessors; ++ii)
	{
		MPI_Status status;
		int nSend = ddc->CommInfo[ii].nSend; 
		int sendStart = ddc->CommInfo[ii].SendStart; 
		int keepSend[nSend+1]; 
	    int remoteId=ddc->CommInfo[ii].domain_id ;
		MPI_Recv(keepSend , nSend , MPI_INT, remoteId, TAG, COMM_LOCAL, &status);
		ddc->CommInfo[ii].SendStart = sendStartNew = ddc->nSend;
		int nSendNew = 0; 
		int *sendlist = ddc->sendList + sendStart;
		int *sendlistNew = ddc->sendList + sendStartNew;
		for (int i = 0; i < nSend; i++) if (keepSend[i]) sendlistNew[nSendNew++] = sendlist[i];
		ddc->CommInfo[ii].SendStart = sendStartNew;
		ddc->CommInfo[ii].nSend = nSendNew;
		ddc->nSend   += nSendNew;
	}
	ddc->number_remote = nremoteNew;
	ddc->number_particles = nremoteNew + nlocal;
	MPI_Status status[nRemoteProcessors];
	MPI_Waitall(nRemoteProcessors, sendRequest, status);
   profile(WAIT_FPAIRS, END);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
