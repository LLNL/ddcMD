#ifndef OMPUTILS_H
#define OMPUTILS_H

#include <omp.h>

#ifdef BGP 
#define CACHELINE_SIZE 128
#else
#define CACHELINE_SIZE 64
#endif


typedef struct cacheline_omp_lock_st
{
   #ifdef WITH_OMP
   omp_lock_t lock;
   char padding[CACHELINE_SIZE-sizeof(omp_lock_t)];
   #else 
   int lock;
   char padding[CACHELINE_SIZE-sizeof(int)];
   #endif

} CACHELINE_OMP_LOCK;

typedef struct cacheline_int_st
{
   double val;
   char padding[CACHELINE_SIZE-sizeof(int)];
} CACHELINE_INT;

typedef struct cacheline_float_st
{
   float val;
   char padding[CACHELINE_SIZE-sizeof(float)];
} CACHELINE_FLOAT;

typedef struct cacheline_double_st
{
   double val;
   char padding[CACHELINE_SIZE-sizeof(double)];
} CACHELINE_DOUBLE;


int getOmpThreadId();
double getOmpWtime();

void setOmpNumThreads(int n);
int  getOmpNumThreads();

void initOmpLock(CACHELINE_OMP_LOCK *lck);
void setOmpLock(CACHELINE_OMP_LOCK *lck);
void unsetOmpLock(CACHELINE_OMP_LOCK *lck);
void destroyOmpLock(CACHELINE_OMP_LOCK *lck);
int testOmpLock(CACHELINE_OMP_LOCK *lck);
#endif

