#ifndef DDC_H
#define DDC_H
#include <mpi.h>
#include <stdio.h>
#include "three_algebra.h"
#include "error.h"
#include "ptiming.h"
#include "external.h"
#include "domain.h"
#include "recordDecorator.h"
#include "gid.h"
#include "loadBalance.h"
#include "ddcRule.h"

#ifdef __cplusplus
extern "C" {
#endif

   //#define MAXPROC 72*32*32*4
#define MAXPROC 96*1024*64
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define NO_DDC_UPDATE   -174321321
#define MSG(text)   location(text,__FILE__,__LINE__)
#ifdef ANSI_F77
#define RINT lrint
#endif
#ifdef AIX
#define LRINT  itrunc
#endif
#ifndef LRINT
#define LRINT lrint
   /*int lrint(double); */
#endif
   char *location(char *, char *, int);
   extern FILE *ddcfile;

   enum DCC_PAIRMODE { ALLREMOTE = 0, LABEL = 1, PROCESSOR = 2 };
   enum DCC_ENUMS { FIRST, NEXT, PARTICLEFUNCTION, BOUNDARY_CONDITION, H, HINV, 
      CORNER, RCUT, DDCNLOCAL, PAIRMODE, GEO, TIMING, CENTER, 
      TWOSIDED, ONESIDED };
   enum DDC_LATTICETYPE { NONE, CUBIC, BCC, FCC, HEXAGONAL, NONLATTICE, RECURSIVE_BISECTION }; 
   enum DDC_TYPE { ARBITRARY, FIXED };

   typedef struct ddcbox_st
   {
      THREE_VECTOR r;
      int nn, first, firstlocal;
      struct ddcbox_st *nlist[256];
   } DDCBOX;

   typedef struct particle_st
   {
      int domain_id;
      int type;
      gid_type global_index;
      unsigned int local_index; 
      unsigned int ifirst;
      THREE_VECTOR r;
   } PARTICLE;


   typedef struct ddctiming_st
   {
      double elapsestart, cpustart, elapsetime, cputime, elapsetotal, cputotal, elapseave, cpuave, ncalls;
      char *name;
   } DDCTIMING;

   typedef struct sendinfo_st
   {
      int domain_id; // task to which data is sent and/or received
      int nSend, SendStart;
      int nRecv, RecvStart;
      void *RecvPtr;
      MPI_Datatype RecvType;
      MPI_Datatype ForceType;
   } COMMINFO;

   typedef struct ddc_st
   {
      char *name;		/* name of the system */
      char *objclass;
      char *value;
      char *type;		/* type label */
      void  *parent; 
      int itype; 
      THREE_VECTOR *corner;
      PARTICLE *particles;
      DOMAINSET domains;
      COMMINFO *CommInfo;  // list of COMMINFO for my neighbors
      
      DDCRULE* ddcRule; // Rule for refine the domain decomposition
      //PTIMING *timing;
      int Assignment_Called; 
      int PairMode;
      int *sendList;
      int domain_id;
      int size;
      enum DDC_LATTICETYPE lattice; // lattice type for subdomains centers
      //	int l;
      int lx,ly,lz;
      double *dx,*dy,*dz; 
      unsigned  number_local, number_remote, number_particles;
      int comm_mode ; 
      double min_NonOverLappingR,min_OverLappingR,max_OverLappingR; 
      int nneighbors,*neighbors;
      int nOverLap;
      int *OverLapID;
      int *startRemoteID; 
      int *nRemoteID; 
      int nRemoteProcessors;
      int nSend;
      int update, updateRate, lastUpdate;
      LOAD_BALANCE* loadBalance;
      double rcut;
      double rInfluence;
      int pruneSends;
      int forceUsePlanes;
      int noDomainSetInAssignmentHack;
      unsigned nPxyzDecorator;
      RECORD_DECORATOR** pxyzDecorator;
      unsigned centersAffineUpdateIndex;
      unsigned pinnedSize;
      unsigned pinnedCapacity;
      unsigned* pinned;
   }
   DDC;

   typedef struct list_st
   {
      void *particle;
      int boxindex, nextindex;
   }
   LIST;

   int ddc_put(DDC*ddc, int put, ...);
   DDC *ddc_init(void *parent,char *name);
   // unused
   //double ddc_get_radius(DDC *ddc0) ;
   THREE_VECTOR ddc_get_center(DDC *ddc0) ;
   // /* Appears not to be used */ void ddc_set_center(THREE_VECTOR new_center, DDC *ddc0); 
   int ddcGetParticleDomainByIndex(DDC *ddc, unsigned i); 
   DDC *getddc(void);
   void ddc_setDomainLattice(DDC* ddc);
   void ddc_setUniformDomainCenters(DDC* ddc);
   void ddc_setRandomDomainCenters(DDC* ddc);
   void ddc_sendLocalRadiusToNeighbors(DDC* ddc);
   void ddc_sendLocalCenterToNeighbors(DDC* ddc, THREE_VECTOR* mycenter);
   void ddc_getDomainSize(DDC*ddc, THREE_VECTOR r0);
   void ddc_initDomains(DDC* ddc);
   void ddc_addPxyzDecorator(DDC* ddc, RECORD_DECORATOR* d);
   void ddc_pinParticle(DDC* self, unsigned index);
   int getLocalTaskOfParticle(unsigned iLocal);
   void ddcUpdate(DDC*ddc);
   void ddcUpdateForce(DDC*ddc,unsigned calculate);
   void ddcGetParticles(DDC*ddc);
   void ddcPutParticles(DDC*ddc);
   void ddcAssignment(DDC*ddc);
   void ddcSendRecvTables(DDC*ddc);

#ifdef __cplusplus
}
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
