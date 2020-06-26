#ifndef THREAD_SCHEDULER_H
#define THREAD_SCHEDULER_H

#include "object.h"
#include "opt.h"
#include "ompUtils.h"

typedef struct schedule_st
{
   int **workalloc;
   int *lastalloc;
   int **deplist;
   int *lastdep;
   int *pcount;
   int **pfmap;
   int nThreads;     // needed for destructor
   int nIdleThreads;
} SCHEDULE;

typedef struct thread_scheduler_st
{
   int nThreads;
   int scheduleRate;
   int fullNN;
   int strictSchedule;
   int profileThread; 
   double imbalanceThreshold; 
   //double* threadTime;
   //double* boxtime;
   CACHELINE_DOUBLE* threadTime; 
   CACHELINE_DOUBLE* boxtime; 
   OPT_BOX_SIZE lastBoxSize;
   SCHEDULE* schedule;
   SCHEDULE* (*buildSchedule) (struct thread_scheduler_st*,
			       OPT_BOX*, OPT_BOX_SIZE);
} THREAD_SCHEDULER;


THREAD_SCHEDULER* threadScheduler_init(OBJECT* obj, int nthreads);

int ts_newScheduleNeeded(THREAD_SCHEDULER* self, OPT_BOX_SIZE box_size, SIGNED64 loop);
int ts_buildSchedule(THREAD_SCHEDULER* self, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize);
void ts_recordBoxTime(THREAD_SCHEDULER* self, int boxid, double time);
void ts_recordThreadTimes(THREAD_SCHEDULER* self, double* time); 
double ts_getMeasuredImbalance(THREAD_SCHEDULER* self);
void ts_mapThreadMemory(THREAD_SCHEDULER* ts, OPT_BOX* boxlist,int tid);

int get_pfmap_offset(int* deplist, int target, int size);

#endif
