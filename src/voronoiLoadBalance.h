#ifndef VORONOI_LOADBALANCE_H
#define VORONOI_LOADBALANCE_H

#include <stdio.h>
#include "loadBalance.h"


void* voronoiLoadBalance_init(LOAD_BALANCE* loadBalance);
void voronoiLoadBalanceFunction(LOAD_BALANCE* this, struct system_st* sys, struct ddc_st* ddc);


void voronoiLoadBalanceWrite(LOAD_BALANCE* this, SIGNED64 loop, double stime, FILE*file);

#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
