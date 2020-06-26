#include "ompUtils.h"

#ifdef WITH_OMP
#include <omp.h>
#endif

int getOmpThreadId()
{
#ifdef WITH_OMP
   return omp_get_thread_num();
#else
   return 0;
#endif
}

double getOmpWtime()
{
#ifdef WITH_OMP
   return omp_get_wtime();
#else
   return 0;
#endif
}

void setOmpNumThreads(int n)
{
#ifdef WITH_OMP
   omp_set_num_threads(n);
#endif
   return;
}

int getOmpNumThreads()
{
#ifdef WITH_OMP
   return omp_get_num_threads();
#else
   return 1;
#endif
}

void initOmpLock(CACHELINE_OMP_LOCK *lck)
{
#ifdef WITH_OMP
   omp_init_lock(&lck->lock);
#else
   return;
#endif
}

void setOmpLock(CACHELINE_OMP_LOCK *lck)
{
#ifdef WITH_OMP
   omp_set_lock(&lck->lock);
#else
   return;
#endif
}

void unsetOmpLock(CACHELINE_OMP_LOCK *lck)
{
#ifdef WITH_OMP
   omp_unset_lock(&lck->lock);
#else
   return;
#endif
}

void destroyOmpLock(CACHELINE_OMP_LOCK *lck)
{
#ifdef WITH_OMP
   omp_destroy_lock(&lck->lock);
#else
   return;
#endif
}

int testOmpLock(CACHELINE_OMP_LOCK *lck)
{
#ifdef WITH_OMP
   return omp_test_lock(&lck->lock);
#else
   return 1;
#endif
}

