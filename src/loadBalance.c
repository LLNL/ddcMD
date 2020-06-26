#include "loadBalance.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "zRampLoadBalance.h"
#include "voronoiLoadBalance.h"

#include "bisectionLoadBalance.h"

#include "object.h"
#include "mpiUtils.h"
#include "ddcMalloc.h"
#include "ddc.h"
#include "gid.h"

// We really don't want to include system.h and neighbor.h, but at the
// moment it is the easiest way to get the number of pairs.  Once we
// have a better solution for tracking npairs we can get rid of the
// coupling.  (This problem used to be handled by storing a copy of the
// NBR pointer in the DDC and getting npair through this pointer.
// However this was the only use of NBR in DDC and we want to avoid
// coupling NBR to DDC if possible.
#include "system.h"
#include "neighbor.h"

void defaultBalanceFunction(LOAD_BALANCE* this, struct system_st* sys, struct ddc_st* ddc){return;};
void defaultWriteFunction(LOAD_BALANCE* this, SIGNED64 loop, double time, FILE* file){return;};


LOAD_BALANCE* loadBalance_init(void* parent, char* name)
{
   LOAD_BALANCE* loadBalance = (LOAD_BALANCE*)
      object_initialize(name, "LOADBALANCE", sizeof(LOAD_BALANCE));
   
   loadBalance->parent = parent;
   loadBalance->rate = 0;
   loadBalance->balanceFunction = defaultBalanceFunction;
   loadBalance->writeFunction = defaultWriteFunction;
   loadBalance->parms = NULL;

   OBJECT* obj = (OBJECT*) loadBalance;
   
   char* type;
   object_get(obj, "rate", &loadBalance->rate, INT,    1, "0");
   object_get(obj, "type", &type,              STRING, 1, "voronoi");
   if (getRank(0) == 0) 
   {
      printf("obj=%s\n",obj->value); 
      printf("loadbalance = %s\n",type); 
   }
   
   if(getRank(0) == 0)
     printf("In %s(): type = %s, rate = %d\n",__func__,type,loadBalance->rate);

   if (strcasecmp(type, "zRamp") == 0 )
   {
      loadBalance->itype = ZRAMP; 
      loadBalance->parms           = zRampLoadBalance_init(loadBalance);
      loadBalance->balanceFunction = zRampLoadBalanceFunction;
   }
   else if (strcasecmp(type, "voronoi") == 0 )
   {
      loadBalance->itype = VORONOI; 
      loadBalance->parms           = voronoiLoadBalance_init(loadBalance);
      loadBalance->balanceFunction = voronoiLoadBalanceFunction;
      loadBalance->writeFunction   = voronoiLoadBalanceWrite;
   }
   else if (strcasecmp(type, "bisection") == 0 )
   {
      loadBalance->itype = BISECTION; 
      loadBalance->parms           = bisectionLoadBalance_init(loadBalance);
      loadBalance->balanceFunction = bisectionLoadBalanceFunction;
   }
   else if (getRank(0) == 0)
   {
      printf("Error: Unsupported load balance type %s\n"
	     "  Game Over\n", type);
      abortAll(1);
   }
   ddcFree(type);
   return loadBalance;
}

double getLocalLoadElapsetime(enum LDBL_LOAD loadMetric, DDC* ddc)
{
   double time=0;
   switch( loadMetric )
   {
     case LDBL_NLOCAL:
      time = (double)ddc->number_local;
      break;
      
     case LDBL_NPAIRS:
      {
	 SYSTEM* sys = system_getSystem(NULL);
	 time = (double)(sys->neighbor->npairs);
      }
      break;
      
     case LDBL_BARRIER:
      time = profileGetLocalLoadElapsetime(0);
      break;
      
     case LDBL_PFORCE:
      time = profileGetLocalLoadElapsetime(1);
      break;
      
     case LDBL_DDCENERGY:
      time = profileGetLocalLoadElapsetime(2);
      break;
      
     case LDBL_REPRODUCIBLE_TEST:
      time = profileGetLocalLoadElapsetime(-1);
      break;
      
     default:
      printf("ERROR - unspecified load measure for load balancing\n");
      assert(1==0);
   }
   
   return time;
}

