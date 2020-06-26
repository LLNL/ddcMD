#include "loadBalance.h"
#include "bisectionLoadBalance.h"

#include <assert.h>

#include "ddcMalloc.h"
#include "object.h"
#include "system.h"
#include "ddc.h"

#include "mpiUtils.h"

void * bisectionLoadBalance_init(LOAD_BALANCE *balancer) {
  BISECTION_BALANCE_PARMS *parms =
    (BISECTION_BALANCE_PARMS *) ddcMalloc(sizeof(*parms));
  OBJECT *obj = (OBJECT *) balancer;
  object_get(obj,"particleWeight",&parms->particle_weight,DOUBLE,1,"1.0");
  object_get(obj,"neighborWeight",&parms->neighbor_weight,DOUBLE,1,"0.0");
  object_get(obj,"volumeWeight",&parms->volume_weight,DOUBLE,1,"0.0");

  parms->reassign = 1;
  if(getRank(0) == 0) 
  {
     printf("In %s(): particle_weight = %.3e, neighbor_weight = %.3e volume_weight = %.3e\n",
           __func__,parms->particle_weight,parms->neighbor_weight,parms->volume_weight);
  }
  const double p = parms->particle_weight,n = parms->neighbor_weight;
  assert(p >= 0.0 && n >= 0.0 && (p+n) >= 0.0);
  assert((p+n)>0 || parms->volume_weight > 0);
  return parms;
}

void bisectionLoadBalanceFunction(LOAD_BALANCE *this,SYSTEM *sys,DDC *ddc) {
   BISECTION_BALANCE_PARMS *parms =
      (BISECTION_BALANCE_PARMS *) (this->parms);
   parms->reassign = 1;
}
void bisectionReAssign(LOAD_BALANCE *loadBalance) 
{
   BISECTION_BALANCE_PARMS *parms = (BISECTION_BALANCE_PARMS *)(loadBalance->parms);
   assert(parms != NULL);
   parms->reassign = 1;
}
