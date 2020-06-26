#ifndef LOAD_BALANCE_H
#define LOAD_BALANCE_H

#include <stdio.h>
#include "gid.h"

struct system_st;
struct ddc_st;

enum LOADBALANCE_CLASS { LB_NONE, VORONOI, ZRAMP, BISECTION}; 
enum LDBL_LOAD { LDBL_NLOCAL, LDBL_NPAIRS, LDBL_BARRIER, LDBL_PFORCE, LDBL_DDCENERGY, LDBL_REPRODUCIBLE_TEST };

typedef struct loadBalance_st
{
   char* name;
   char* objclass;
   char* value;
   char* type;
   void* parent; 
   enum LOADBALANCE_CLASS itype; 
   int rate;
   void (*balanceFunction) (struct loadBalance_st* self,
			    struct system_st* sys,
			    struct ddc_st* ddc);
   void (*writeFunction)   (struct loadBalance_st* self,
			    SIGNED64 loop,
			    double time,
			    FILE* file);
   void* parms;
} LOAD_BALANCE;

LOAD_BALANCE* loadBalance_init(void* parent, char* name);

double getLocalLoadElapsetime(enum LDBL_LOAD loadMetric, struct ddc_st* ddc);

#endif
