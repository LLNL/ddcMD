#include "ddc.h"
#include <assert.h>
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
#include "object.h"
#include "particle.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "preduce.h"
#include "readPXYZ.h"
#include "utilities.h"
#include "domain.h"

#define TAG 1

//int domain_id;

static DDC *_ddc = NULL;

void get_particles(PARTICLE*ddc_particle, int  ddc_n);
void put_particles(PARTICLE*ddc_particle, int  nlocal, int  nremote);
int ddc_put(DDC*ddc, int put, ...);
void domain_set(DDC*ddc);



DDC *ddc_init(void *parent,char *name)
{
   DDC *ddc;
   int i,size,nelements;
   char *string=NULL; 
   char *type=NULL; 
   double d[1024]; 
   ddc = (DDC *) object_initialize(name, "DDC", sizeof(DDC));
   ddc->parent = parent; 

   MPI_Comm_size(COMM_LOCAL, &size);
   ddc->size = size;
   MPI_Comm_rank(COMM_LOCAL, &ddc->domain_id);

   nelements = object_get((OBJECT *) ddc, "type", &type, STRING, 1, "ARBITRARY");
   ddc->type = ddcCalloc(strlen(type) + 1, sizeof(char)); strcpy(ddc->type, type);
   if (strcmp(type, "ARBITRARY") == 0) ddc->itype = ARBITRARY;
   if (strcmp(type, "FIXED"    ) == 0) ddc->itype = FIXED;

   int lx,ly,lz;
   nelements = object_get((OBJECT *) ddc, "lx", &lx, INT, 1, "0");
   nelements = object_get((OBJECT *) ddc, "ly", &ly, INT, 1, "0");
   nelements = object_get((OBJECT *) ddc, "lz", &lz, INT, 1, "0");

   ddc->loadBalance = NULL;
   char* balName = NULL;
   object_get((OBJECT *) ddc, "loadbalance", &(balName), STRING, 1, "loadbalance");
   if (getRank(0)==0) 
   {
      printf("ddc=%s\n",ddc->value); 
      printf("balName=%s\n",balName); 
   }
   ddc->loadBalance = loadBalance_init(ddc, balName);
   ddcFree(balName);

   // read relative domain centers
   ddc->dx=ddc->dy=ddc->dz=NULL; 
   nelements = object_get((OBJECT *)ddc,"dx",d,DOUBLE,1024,"-1.0"); 
   if (d[0] > 0.0)  { lx = nelements ;  ddc->dx = ddcMalloc(lx*sizeof(double));for (i=0;i<lx;i++) ddc->dx[i] = d[i];}

   nelements = object_get((OBJECT *)ddc,"dy",d,DOUBLE,1024,"-1.0"); 
   if (d[0] > 0.0)  { ly = nelements ;  ddc->dy = ddcMalloc(ly*sizeof(double));for (i=0;i<ly;i++) ddc->dy[i] = d[i];}

   nelements = object_get((OBJECT *)ddc,"dz",d,DOUBLE,1024,"-1.0"); 
   if (d[0] > 0.0)  { lz = nelements ;  ddc->dz = ddcMalloc(lz*sizeof(double));for (i=0;i<lz;i++) ddc->dz[i] = d[i];}
   if (ddc->size == 1 ) lx = ly = lz=1; 

   ddc->lx=lx; 
   ddc->ly=ly; 
   ddc->lz=lz; 

   // get initial lattice type for domain decomposition
   nelements = object_get((OBJECT *) ddc, "lattice", &string, STRING, 1, "NONE"); 
   ddc->lattice = NONE; 
   if (!strcmp(string,"CUBIC")) ddc->lattice = CUBIC; 
   if (!strcmp(string,"BCC")) ddc->lattice = BCC; 
   if (!strcmp(string,"FCC")) ddc->lattice = FCC; 
   if (!strcmp(string,"HEXAGONAL")) ddc->lattice = HEXAGONAL; 
   if (!strcmp(string,"NONLATTICE")) ddc->lattice = NONLATTICE; 
   if (!strcmp(string,"RECURSIVE_BISECTION")) ddc->lattice = RECURSIVE_BISECTION;

   nelements = object_get((OBJECT *) ddc, "updateRate", &ddc->updateRate, INT, 1, "0");
   nelements = object_get((OBJECT *) ddc, "pruneSends", &ddc->pruneSends, INT, 1, "1");
   nelements = object_get((OBJECT *) ddc, "forceUsePlanes", &ddc->forceUsePlanes, INT, 1, "0");
   object_get((OBJECT *) ddc, "noDomainSetInAssignmentHack", &ddc->noDomainSetInAssignmentHack, INT, 1, "0");
   if (size > MAXPROC)
   {
      printf("Exceeded number of allowed Processes. Increase MAXPROC in ddc.h and Recompile\n"),
         abortAll(-1); 
   }
   
   char *ddcRuleName;
   nelements = object_get((OBJECT *) ddc, "ddcRule", &ddcRuleName, STRING, 1, "defaultRue");
   if(strcmp(ddcRuleName, "defaultRue")==0){
       ddc->ddcRule=NULL;
   }else{
       ddc->ddcRule=ddcRule_init(ddc, ddcRuleName);
   }
   ddcFree(ddcRuleName);

   ddc->Assignment_Called = 0;
   ddc->number_particles = 0;
   ddc->particles = NULL;
   ddc->sendList = NULL;
   ddc->CommInfo = NULL;
   ddc->nRemoteProcessors = 0;
   ddc->update=1;
   ddc->lastUpdate=NO_DDC_UPDATE;
   ddc->nOverLap=0; 
   ddc->OverLapID=NULL;
   ddc->nPxyzDecorator = 0;
   ddc->pxyzDecorator = NULL;

   ddc->pinnedSize = 0;
   ddc->pinnedCapacity = 0;
   ddc->pinned = NULL;

   ddc_initDomains(ddc);

   _ddc = ddc;

   return ddc;
}
void ddc_initDomains(DDC* ddc)
{
   timestamp("Start ddc_initDomains");

   ddc_setDomainLattice( ddc );
   domainset_init( &ddc->domains );

   int tryReadPxyz = (ddc->lattice == NONLATTICE);
   // If you're hitting this assertion I owe you an apology.  Sorry.  On
   // 29-Aug-2011 it was a high priority to get ddc.c decoupled from
   // SIMULATE.  This was the last remaining coupling.  It was more
   // important to get the decoupling done that to keep the pxyz reading
   // working for load balancing and I didn't have time to do both.
   //
   // I don't remember why it was important that we not attempt to read
   // pxyz files if the load balance rate was non-zero.  However,
   // without the SIMULATE object I can't make the check.  To get up and
   // running again it will be necessary to revisit the logic for
   // reading pxyz files and establish once and for all what the proper
   // behavior is.  Whatever we decide, we have to do it without adding
   // a requirement to access simulate.
   //
   // A little more information (13-Dec-2011)
   // I just fixed a bug in r1444 for Stephane Mazevet.  In that
   // version, reading pxyz unconditionally sets ddc->lattice to
   // NONLATTICE.  This effectively disables the zRampLoadBalance, since
   // it relies on the ability to call domain_set and get lattice-like
   // behavior.  Note that this version of the code won't even try to
   // read pxyz in the case that you're doing zRamp.  That probably
   // isn't what we want either.  I'm going to need to think hard about this.

   //ewd: readPxyz is useful for Voronoi runs on foams.  We're not currently using
   //ewd:  zRampLoadBalance, so I'm disabling this assert for now
   //ewd	if (tryReadPxyz == 1) assert(1==0);
   /* 	{  // limit scope of SIMULATE.  We shouldn't depend on SIMULATE. */
   /* 		SIMULATE* simulate = simulate_getSimulate(NULL); */
   /* 		tryReadPxyz |= (simulate->loadBalance->rate != 0); */
   /* 	} */
   int readStatus = -1;
   if ( tryReadPxyz )
      readStatus = readPXYZ(ddc);

   if ( readStatus == 0 )
      ddc->lattice = NONLATTICE;
   else // read failed (or wasn't attempted)
   {
      if ( ddc->lattice==NONLATTICE )
         domainset_randomCenters(&ddc->domains, ddc->lx, ddc->ly, ddc->lz);
      else
         domain_set(ddc);
   }
   timestamp("End ddc_initDomains");
}
void ddcPutParticles(DDC*ddc)
{
   put_particles(ddc->particles, ddc->number_local, ddc->number_remote);
}
int ddcGetParticleDomainByIndex(DDC *ddc,unsigned index)
{
   unsigned maxIndex=ddc->number_local;
   if (index < maxIndex) return ddc->domain_id;

   for (int rIndex =0;rIndex<ddc->nRemoteProcessors;rIndex++)
   {
      maxIndex += ddc->CommInfo[rIndex].nRecv;
      if (index < maxIndex) return ddc->CommInfo[rIndex].domain_id; 
   }
   return -1; 
}

void ddcGetParticles(DDC*ddc)
{
   get_particles(ddc->particles, ddc->number_local);
}

DDC *getddc(void)
{
   return  _ddc; 
}

//unused
/* double ddc_get_radius(DDC *ddc0)  */
/* { */
/*  	if (ddc0==NULL)  ddc0 = _ddc;  */
/*    return domainset_getLocalRadius( &ddc0->domains ); */
/* } */
THREE_VECTOR ddc_get_center(DDC *ddc0) 
{
   if (ddc0==NULL)  ddc0 = _ddc; 
   return domainset_getLocalCenter( &ddc0->domains ); 
}
/* This is never used. Let's comment it out, since this call
   probably requires a following call to domainset_allGather() */
/*
   void ddc_set_center(THREE_VECTOR new_center, DDC *ddc0) 
   {
   if (ddc0==NULL)  ddc0 = _ddc; 
   domainset_setLocalCenter( &ddc0->domains, new_center);
   }
   */

int ddc_put(DDC*ddc, int put, ...)
{
   va_list ap;
   va_start(ap, put);

   switch (put)
   {
      case CORNER:
         ddc->corner = va_arg(ap,THREE_VECTOR *);
         va_end(ap); 
         return 0;
      case DDCNLOCAL:
         /*ddc->number_local = (int)ptr;*/
         ddc->number_local = va_arg(ap,unsigned);
         va_end(ap); 
         return 0;
      case PAIRMODE:
         ddc->PairMode = va_arg(ap,int);
         va_end(ap); 
         return 0;
      case RCUT:
         ddc->rcut = va_arg(ap,double );
         va_end(ap); 
         return 0;
      default:
         va_end(ap); 
         break;
   }
   return 1;
}

// unused function
/* double ddcweightpair_(int *iptr, int *jptr) */
/* { */
/* 	unsigned int ii, ij;  */
/* 	int i, j; */
/* 	i = (*iptr) - 1; */
/* 	j = (*jptr) - 1; */
/* 	ii = _ddc->particles[i].global_index; */
/* 	ij = _ddc->particles[j].global_index; */
/* 	if (j < _ddc->number_local) */
/* 	{ */
/* 		if (ij > ii) return 1.0; */
/* 		return 0; */
/* 	} */
/* 	switch (_ddc->PairMode) */
/* 	{ */
/* 	case ALLREMOTE: */
/* 		return 0.5; */
/* 	case LABEL: */
/* 		if (ij > ii) return 1.0; */
/* 		return 0.0; */
/* 	case PROCESSOR: */
/* 		if (_ddc->particles[j].domain_id > _ddc->particles[i].domain_id) return 1.0; */
/* 		return 0.0; */
/* 	} */
/* 	return -1.0;  */
/* } */

void ddc_setUniformRelativeDomainCenters(DDC* ddc)
{
   if ( ddc->dx == NULL || ddc->dy == NULL || ddc->dz == NULL ){
      if (ddc->domain_id == 0)
         printf("ddc_setUniformRelativeDomainCenters for lx=%d, ly=%d, lz=%d...\n",
               ddc->lx, ddc->ly, ddc->lz);
   }
   if ( ddc->dx == NULL ) {
      ddc->dx = ddcMalloc(ddc->lx*sizeof(double));
      for (int i=0;i<ddc->lx;i++) ddc->dx[i] = (1.0/ddc->lx )*(i + 0.5);
   }
   if ( ddc->dy == NULL ) {
      ddc->dy = ddcMalloc(ddc->ly*sizeof(double));
      for (int i=0;i<ddc->ly;i++) ddc->dy[i] = (1.0/ddc->ly )*(i + 0.5);
   }
   if ( ddc->dz == NULL ) {
      ddc->dz = ddcMalloc(ddc->lz*sizeof(double));
      for (int i=0;i<ddc->lz;i++) ddc->dz[i] = (1.0/ddc->lz )*(i + 0.5);
   }
}

void ddc_setDomainLattice(DDC* ddc)
{
   timestamp("Start ddc_setDomainLattice");
   if( ddc->lattice != NONLATTICE  &&  ddc->lattice != RECURSIVE_BISECTION ){
      int n = ddc->size; // # MPI tasks
      long lmax = LRINT(ceil(cbrt(1.0*MAXPROC)));
      // determine regular lattice type

      if (ddc->lx*ddc->ly*ddc->lz==0 )
      {
         int lattice_type = NONE;
         int l;
         for (l = 0; l < lmax; l++)
         {
            if (ddc->lattice == NONE) 
            {
               if (  l*l*l == n) lattice_type = CUBIC;
               if (2*l*l*l == n) lattice_type = BCC;
               if (4*l*l*l == n) lattice_type = FCC;
            }
            else
            {
               if (  l*l*l == n && ddc->lattice == CUBIC) lattice_type = CUBIC;
               if (2*l*l*l == n && ddc->lattice == BCC)   lattice_type = BCC;
               if (4*l*l*l == n && ddc->lattice == FCC)   lattice_type = FCC;
               if (2*l*l*l == n && ddc->lattice == HEXAGONAL) lattice_type = HEXAGONAL;
            }
            if (lattice_type != NONE) break;
         }
         ddc->lx=ddc->ly=ddc->lz = l;
         ddc->lattice = lattice_type;
      }
      else 
      {
         int lx = ddc->lx;
         int ly = ddc->ly;
         int lz = ddc->lz;
         if (ddc->lattice == NONE  )
         {
            int lattice_type = NONE;
            if (  lx*ly*lz == n) lattice_type = CUBIC;
            if (2*lx*ly*lz == n) lattice_type = BCC;
            if (4*lx*ly*lz == n) lattice_type = FCC;
            ddc->lattice = lattice_type;
         }

      }

      ddc_setUniformRelativeDomainCenters(ddc);
   }
}

void ddc_sendLocalRadiusToNeighbors(DDC* ddc)
{
   assert( ddc->nRemoteProcessors>0 );

   MPI_Request* recvRequest = ddcMalloc(ddc->nRemoteProcessors*sizeof(MPI_Request));
   MPI_Request* sendRequest = ddcMalloc(ddc->nRemoteProcessors*sizeof(MPI_Request));

   double send_buffer[1]={domainset_getLocalRadius(&ddc->domains)};
   double* recv_buffer=ddcMalloc(ddc->nRemoteProcessors*sizeof(double));
   const int tag=19;
   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      assert( remoteID!=ddc->domain_id );
      MPI_Irecv(recv_buffer+i, 1, MPI_DOUBLE, remoteID, tag, COMM_LOCAL, recvRequest + i);
   }

   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      MPI_Isend(send_buffer, 1, MPI_DOUBLE, remoteID, tag, COMM_LOCAL, sendRequest + i);
   }

   MPI_Waitall(ddc->nRemoteProcessors, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(ddc->nRemoteProcessors, recvRequest, MPI_STATUSES_IGNORE);

   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      domainset_setRadius( &ddc->domains, remoteID, recv_buffer[i]);
   }

   ddcFree(recv_buffer);

   ddcFree(recvRequest);
   ddcFree(sendRequest);
}
void ddc_sendLocalCenterToNeighbors(DDC* ddc, THREE_VECTOR* mycenter)
{
   MPI_Request* recvRequest = ddcMalloc(ddc->nRemoteProcessors*sizeof(MPI_Request));
   MPI_Request* sendRequest = ddcMalloc(ddc->nRemoteProcessors*sizeof(MPI_Request));

   double send_buffer[3]={mycenter->x,mycenter->y,mycenter->z};
   double* recv_buffer=ddcMalloc(ddc->nRemoteProcessors*3*sizeof(double));
   const int tag=18;
   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      assert( remoteID!=ddc->domain_id );
      MPI_Irecv(recv_buffer+3*i, 3, MPI_DOUBLE, remoteID, tag, COMM_LOCAL, recvRequest + i);
   }  

   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      MPI_Isend(send_buffer, 3, MPI_DOUBLE, remoteID, tag, COMM_LOCAL, sendRequest + i);
   }

   MPI_Waitall(ddc->nRemoteProcessors, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(ddc->nRemoteProcessors, recvRequest, MPI_STATUSES_IGNORE);

   for (int i = 0; i < ddc->nRemoteProcessors; i++)
   {
      int remoteID = ddc->CommInfo[i].domain_id;
      VSET(ddc->domains.domains[remoteID].center,
            recv_buffer[3*i],recv_buffer[3*i+1],recv_buffer[3*i+2]);
   }  
   ddcFree(recv_buffer);

   ddcFree(recvRequest);
   ddcFree(sendRequest);
}

void ddc_getDomainSize(DDC*ddc, THREE_VECTOR r0)
{
   THREE_VECTOR r;
   THREE_VECTOR extent = vzero;
   double r2max = 0.0;
   for (unsigned i = 0; i < ddc->number_local; i++)
   {
      r = ddc->particles[i].r;
      r.x -= r0.x;
      r.y -= r0.y;
      r.z -= r0.z;
      nearestImage(&r.x, &r.y, &r.z);
      double r2 = (r.x*r.x + r.y*r.y + r.z*r.z);
      if (r2 > r2max) r2max = r2;
      extent.x = MAX(extent.x, fabs(r.x));
      extent.y = MAX(extent.y, fabs(r.y));
      extent.z = MAX(extent.z, fabs(r.z));
   }

   double radius = sqrt(r2max); 
   double volSphere = 4*M_PI/3.0 * CUBE(radius); 
   double volBrick  = extent.x*extent.y*extent.z; 
   int shape = SPHERE ; 
   if (volBrick < volSphere) shape=BRICK; 
   domainset_setLocalShape( &ddc->domains, shape );
   domainset_setLocalRadius( &ddc->domains, radius );
   domainset_setLocalExtent(&(ddc->domains),extent);
}

void ddc_addPxyzDecorator(DDC* this, RECORD_DECORATOR* d)
{
   this->pxyzDecorator = ddcRealloc(
         this->pxyzDecorator, (1+this->nPxyzDecorator)*sizeof(RECORD_DECORATOR*));
   this->pxyzDecorator[this->nPxyzDecorator] = d;
   ++this->nPxyzDecorator;
}

void ddc_pinParticle(DDC* this, unsigned index)
{
   if (this->pinnedSize == this->pinnedCapacity)
      this->pinnedCapacity = MAX(10, 2*this->pinnedCapacity);
   this->pinned = ddcRealloc(this->pinned, this->pinnedCapacity*sizeof(unsigned));
   this->pinned[this->pinnedSize++] = index;
}
int getLocalTaskOfParticle(unsigned iLocal)
{
   DDC* ddc = getddc();
   int domain_id = ddc->domain_id; 
   unsigned recvStart=ddc->number_local;
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
   {
      recvStart = ddc->CommInfo[ii].RecvStart ; 
      //int recvEnd = recvStart + ddc->CommInfo[ii].nRecv;
      if   (iLocal <recvStart) return domain_id; 
      domain_id = ddc->CommInfo[ii].domain_id;
   }
   return domain_id; 
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
