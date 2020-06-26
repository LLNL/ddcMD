#ifndef ZRAMP_LOAD_BALANCE_H
#define ZRAMP_LOAD_BALANCE_H

#include "loadBalance.h"

void* zRampLoadBalance_init(LOAD_BALANCE* loadBalance);
void zRampLoadBalanceFunction(LOAD_BALANCE* this, struct system_st* sys, struct ddc_st* ddc);

#endif
