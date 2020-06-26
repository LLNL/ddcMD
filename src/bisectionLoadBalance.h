#ifndef BISECTION_LOADBALANCE__
#define BISECTION_LOADBALANCE__

#include "loadBalance.h"
#include "system.h"
#include "ddc.h"

typedef struct bisectionBalanceParms_st {
  double particle_weight,neighbor_weight,volume_weight;
  int reassign;
} BISECTION_BALANCE_PARMS;

void * bisectionLoadBalance_init(LOAD_BALANCE *balancer);
void bisectionLoadBalanceFunction(LOAD_BALANCE *this,SYSTEM *sys,DDC *ddc);
void bisectionReAssign(LOAD_BALANCE *loadBalance) ;

#endif
