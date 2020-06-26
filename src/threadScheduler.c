#include "threadScheduler.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "utilities.h"
#include "mpiUtils.h"
#include "object.h"
#include "ptiming.h"
#include "ddcMalloc.h"
#include "format.h"

static SCHEDULE* schedule_init(OPT_BOX_SIZE boxsize, int nThreads);
static void schedule_destroy(SCHEDULE* this);

static int ts_optBoxChange(THREAD_SCHEDULER* this, OPT_BOX_SIZE boxSize);
static int ts_improvedSchedule(THREAD_SCHEDULER* this, SCHEDULE* newSchedule);
static void ts_installSchedule(THREAD_SCHEDULER* this, OPT_BOX_SIZE boxSize, SCHEDULE* schedule);
static double ts_maxLoad(THREAD_SCHEDULER* this, SCHEDULE* schedule);
static void ts_estimateBoxTimes(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE box_size);

static void ts_addWork(OPT_BOX *boxlist, int **growthlist, int *lastgrowth, int targetBoxId, int targetThread, OPT_BOX_SIZE boxSize, int fullNN, int profileThread); 
static void ts_removeWork(SCHEDULE *sc, OPT_BOX *boxlist, int **growthlist, int *lastgrowth, double *threadWorkload, double targetBoxtime, int targetBoxId, int targetThread, OPT_BOX_SIZE boxSize, int fullNN); 
static void ts_addGrowthList(int *growthlist,int *nntable,int *lastdept, int nntablecount);
static void ts_buildThreadGrowthList(SCHEDULE *sc, OPT_BOX *boxlist, int *growthlist, int *lastgrowth, int targetThreadId,OPT_BOX_SIZE boxSize, int fullNN);
static void ts_buildDeptList(SCHEDULE *sc, OPT_BOX *boxlist, int actual_nThreads, int nbox);	

SCHEDULE* buildSchedule_MCAS(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize);
SCHEDULE* buildSchedule_BFAS(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize);
SCHEDULE* buildSchedule_trivial(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize);

static FILE *profileFile;
static SIGNED64 ts_loop;
 
THREAD_SCHEDULER* threadScheduler_init(OBJECT* obj, int nThreads)
{
   THREAD_SCHEDULER* this = ddcMalloc(sizeof(THREAD_SCHEDULER));
   this->boxtime = NULL;
   //this->threadTime = ddcMalloc(nThreads*sizeof(double));
   this->threadTime = ddcMalloc(nThreads*sizeof(CACHELINE_DOUBLE));

   for (int i=0;i<nThreads;i++)
      this->threadTime[i].val = 1.0;

   this->lastBoxSize.nx = 0;
   this->lastBoxSize.ny = 0;
   this->lastBoxSize.nz = 0;
   this->lastBoxSize.nbox = 0;

   this->nThreads = nThreads; 
   this->schedule = NULL;
   this->buildSchedule = NULL;

   char* method;
   object_get(obj, "scheduleRate",       &this->scheduleRate, INT, 1, "0");
   object_get(obj, "imbalanceThreshold", &this->imbalanceThreshold, DOUBLE , 1, "0.0");
   object_get(obj, "scheduleMethod",     &method, STRING, 1, "MCAS");
   object_get(obj, "scheduleFullNN",     &this->fullNN, INT, 1, "1");
   object_get(obj, "strictSchedule",     &this->strictSchedule, INT, 1, "0");
   object_get(obj, "profileThread",     &this->profileThread, INT, 1, "0");

   if (strcasecmp(method, "MCAS")==0) 
      this->buildSchedule=buildSchedule_MCAS;
   if (strcasecmp(method, "BFAS")==0) 
      this->buildSchedule=buildSchedule_BFAS;
   if (nThreads==1) 
      this->buildSchedule=buildSchedule_trivial;

   


   char combinedName[40];
   char strict[10] = "";
   char fullNN[10] = "";
   if (this->fullNN==1)
      strcpy(fullNN,"-FullNN");
   if (this->strictSchedule==1)
      strcpy(strict,"-Strict");
   
   sprintf(combinedName,"%s%s%s",method,fullNN,strict);

   profileThreadScheduleMethodName(combinedName);

   if (this->strictSchedule == 1 && this->fullNN == 0)
   {
      if (getRank(0)==0) 
         printf("Error: only scheduleFullNN = 1 allowed for strictSchedule mode\n");
      abortAll(1);
   }

   if (this->buildSchedule == NULL) 
   {
      if (getRank(0)==0) 
         printf("Error: Invalid buildSchedule method = %s\n",method);
      abortAll(1);
   }

   return this;
}

/**  There are three reasons we might need a new schedule:
 *   - The size of the OPT_BOX has changed from the previous step (this
 *     invalidates the previous schedule).
 *   - The loop count is a multiple of the schedule rate.
 *   - The imbalance threshhold is exceeded.
 */
int ts_newScheduleNeeded(THREAD_SCHEDULER* this, OPT_BOX_SIZE box_size, SIGNED64 loop)
{
   ts_loop = loop;
   if (ts_optBoxChange(this, box_size))
      return 1;
   if (TEST0(loop, this->scheduleRate))
      return 1;
   if ( (this->imbalanceThreshold > 0.0) &&
        (ts_getMeasuredImbalance(this) > this->imbalanceThreshold) )
      return 1;
   return 0;
}

/** This function doesn't quite do what you might think since you may
 *  not get a new schedule when you call it.
 *
 *  It does always build a new schedule, but it only "installs" it in
 *  the THREAD_SCHEDULER if it is "better" than the schedule that is
 *  already there.
 *
 *  The new schedule is better if the OPT_BOX_SIZE size has changed
 *  (this makes the old schedule invalid) or if the ts_improvedSchedule
 *  function says so.
 */
int ts_buildSchedule(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize)
{
   int boxChange = ts_optBoxChange(this, boxSize);

   if (boxChange)
      ts_estimateBoxTimes(this, boxlist, boxSize);
 
   if (this->profileThread == 1)
   {
      char fname[30];
      sprintf(fname,"%8.8d-",getRank(0));
      sprintf(fname+strlen(fname),loopFormat(),ts_loop);
      sprintf(fname+strlen(fname),".pack");
      profileFile = fopen(fname,"w");
      fprintf(profileFile,"%d\n",this->nThreads);
      fprintf(profileFile,"%d %d %d\n",boxSize.nx,boxSize.ny,boxSize.nz);
   }

  
   profileMT(-4);
   SCHEDULE* newSchedule = this->buildSchedule(this, boxlist, boxSize);

   if (this->profileThread == 1)
   {
      fclose(profileFile);
   }

   if (boxChange || ts_improvedSchedule(this, newSchedule))
      ts_installSchedule(this, boxSize, newSchedule);
   else
      schedule_destroy(newSchedule);
 
   profileNBox((boxSize.nx-2)*(boxSize.ny-2)*(boxSize.nz-2));
   return this->schedule->nIdleThreads;
}

void ts_recordBoxTime(THREAD_SCHEDULER* this, int boxid, double time)
{
   this->boxtime[boxid].val = time;
}

void ts_recordThreadTimes(THREAD_SCHEDULER* this, double* time)
{
   for (int ii=0; ii<this->nThreads; ++ii) 
   {
      this->threadTime[ii].val=time[ii];
   }  
}


double ts_getMeasuredImbalance(THREAD_SCHEDULER* this)
{
   double sumtime=this->threadTime[0].val;
   double maxtime=this->threadTime[0].val;

   for (int i=1; i<this->nThreads; ++i)
   {
      sumtime += this->threadTime[i].val;
      if (maxtime < this->threadTime[i].val) 
         maxtime = this->threadTime[i].val;
   }
   
   return (maxtime/(sumtime/(this->nThreads)))-1.0;
}

//binary search
int get_pfmap_offset(int* deplist, int target,int size)
{
   int mid = size/2;
   int low = 0;
   int high = size-1;
   while (deplist[mid] != target)
   {	
      if (low > high)
      { //this mean search is not found: not allow, all dependency should be mapped 
         #pragma execution_frequency(very_low)
	 printf("Error: Request for non exist element: target=%5d deplist_size=%5d\n",target,size);
	 printf("Deplist elements:");
	 for (int i = 0;i < size;i++)
	    printf("%5d ",deplist[i]);
	 printf("\n");	
	 abortAll(1);
      }	 
      if (deplist[mid] > target)
	 high = mid-1;
      else
	 low = mid+1;
      mid = (low+high)/2;	
   }
   return mid;
}


///// PRIVATE FUNCTIONS /////

SCHEDULE* schedule_init(OPT_BOX_SIZE boxsize, int nThreads)
{
   SCHEDULE* this = ddcMalloc(sizeof(SCHEDULE));
   this->nThreads = nThreads; 
   
   this->workalloc = ddcMalloc(nThreads*sizeof(int *));
   this->lastalloc = ddcMalloc(nThreads*sizeof(int));
   this->deplist = ddcMalloc(nThreads*sizeof(int *));
   this->lastdep = ddcMalloc(nThreads*sizeof(int));
   this->pfmap = ddcMalloc(nThreads*sizeof(int*));
   this->pcount = ddcMalloc(nThreads*sizeof(int));

   for (int i=0;i<nThreads;i++)
   {
      this->lastalloc[i] = 0;
      this->lastdep[i] = 0;
      this->pcount[i] = 0;
      this->workalloc[i] = ddcMalloc(sizeof(int)*boxsize.nbox);
      this->deplist[i] = ddcMalloc(sizeof(int)*boxsize.nbox);
      this->pfmap[i] = NULL;
   }
   return this;
}

void schedule_destroy(SCHEDULE* this)
{
   if (this == NULL) return;  
   for (int i=0; i<this->nThreads; ++i)
   {
      ddcFree(this->workalloc[i]);
      ddcFree(this->deplist[i]);
      ddcFree(this->pfmap[i]);
   }

   ddcFree(this->workalloc);
   ddcFree(this->lastalloc);
   ddcFree(this->deplist);
   ddcFree(this->lastdep);
   ddcFree(this->pfmap);
   ddcFree(this->pcount);
   ddcFree(this);
}

int ts_optBoxChange(THREAD_SCHEDULER* this, OPT_BOX_SIZE boxSize)
{
   if (this->lastBoxSize.nx != boxSize.nx ||
       this->lastBoxSize.ny != boxSize.ny ||
       this->lastBoxSize.nz != boxSize.nz )
      return 1;
   return 0;
}

/** This function checks to see if the new schedule is an improvement
 *  over the schedule currently installed in the the THREAD_SCHEDULER.
 *  The new schedule is considered an improvement if any of the
 *  following are true:
 *
 *  * It has fewer idle threads
 *  * The workload of the most loaded thread is lower
 *
 *  Returns 1 if the new schedule is an improvement, zero otherwise.
 */
int ts_improvedSchedule(THREAD_SCHEDULER* this, SCHEDULE* newSchedule)
{
   if (this->schedule==NULL)
      return 1;
   if (newSchedule->nIdleThreads < this->schedule->nIdleThreads)
      return 1;
   if (ts_maxLoad(this, newSchedule) < ts_maxLoad(this, this->schedule))
      return 1;
   
   return 0;
}

void ts_installSchedule(THREAD_SCHEDULER* this, OPT_BOX_SIZE boxSize, SCHEDULE* schedule)
{
   schedule_destroy(this->schedule);
   this->schedule = schedule;
   this->lastBoxSize = boxSize;
}

/** Computes the load of the most loaded thread of the given schedule,
 *  based on the boxtimes in the given THREAD_SCHEDULER.
 *
 *  Note that this function does not operate on the scheduler stored in
 *  this.*/
double ts_maxLoad(THREAD_SCHEDULER* this, SCHEDULE* schedule) 
{
   if (schedule->nThreads == 1)
      return 0.0;

   double load[schedule->nThreads];

   for (int ii=1; ii<schedule->nThreads; ++ii)
   {
      load[ii] = 0.0;
      for (int jj=0; jj<schedule->lastalloc[ii]; ++jj)
	 load[ii] += this->boxtime[jj].val;
   }

   double maxLoad = load[0];
   for (int ii=1; ii<schedule->nThreads; ++ii)
      if (maxLoad < load[ii]) maxLoad = load[ii];

   return maxLoad;
}

/** Populates the this->boxtime array with estimated workload for each box.
 *  The units of workload are arbitrary as they are only compared with
 *  each other.  */
void ts_estimateBoxTimes(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE box_size)
{
   int nx = box_size.nx;
   int ny = box_size.ny;
   int nz = box_size.nz;
   int nbox = box_size.nbox;
	
   //this->boxtime = ddcRealloc(this->boxtime, nbox*sizeof(double));
   this->boxtime = ddcRealloc(this->boxtime, nbox*sizeof(CACHELINE_DOUBLE));
   for (int iz=1; iz<nz-1; iz++)
   {
      int i =  1 + nx*(1 + ny *iz);
      if (boxlist[i].edge == 0 )
      {
	 for (int iy=1; iy<ny-1; iy++)
	 {
	    for (int ix=1; ix<nx-1; ix++)
	    {
	       int workcount = 0;
	       OPT_BOX *box_i = boxlist+i;
	       int nTotal_i = box_i->nlocal+box_i->nremote;
	       if (nTotal_i <= 0)
	       {
		  i++;
		  continue;
	       }
	       workcount += nTotal_i*(nTotal_i-1);
	       for  (int k=0; k<box_i->nn; k++) /* loop over neighboring boxes of box */
	       {
		  OPT_BOX *box_j = box_i+box_i->nlist[k].offset;
		  int nTotal_j = box_j->nlocal+box_j->nremote;
		  workcount += nTotal_i*nTotal_j;
	       }
	       this->boxtime[i].val=workcount;
	       i++;
	    }
	    i+=2;
	 }
      }
   }
}
/*
void ts_mapThreadMemory(THREAD_SCHEDULER *this, OPT_BOX* boxlist)
{
  SCHEDULE *schedule = this->schedule;
  //counting private allocation space for each thread from dependency list   
  for (int k=0;k < schedule->nThreads;k++)
  {
    schedule->pcount[k] = 0;
    for (int i = 0;i < schedule->lastdep[k];i++)
    {
      OPT_BOX *box_i = boxlist+schedule->deplist[k][i];
      schedule->pcount[k]+= box_i->nlocal+box_i->nremote;
    }
  }

  //define mapping between main force array and truncated version of threadprivate
  for (int k = 0;k < schedule->nThreads;k++)
  {
    schedule->pfmap[k] = ddcRealloc(schedule->pfmap[k], schedule->lastdep[k]*sizeof(int));
    int count = 0;
    for (int i = 0;i < schedule->lastdep[k];i++)
    {
      OPT_BOX *box_tmp = boxlist + schedule->deplist[k][i];
      schedule->pfmap[k][i] = count;
      count+=box_tmp->nlocal+box_tmp->nremote;
    }
  }
}
*/

//mapThreadMemory by individual thread. Currently not working due to ddcRealloc cannot be called
//by the thread.  
void ts_mapThreadMemory(THREAD_SCHEDULER *this, OPT_BOX* boxlist,int tid)
{
  SCHEDULE *schedule = this->schedule;
  //counting private allocation space for each thread from dependency list   
  //printf("myid=%5d tid = %d schedule=%p\n",getRank(0),tid,schedule);  
  schedule->pcount[tid] = 0;
  for (int i = 0;i < schedule->lastdep[tid];i++)
  {
    OPT_BOX *box_i = boxlist+schedule->deplist[tid][i];
    schedule->pcount[tid]+= box_i->nlocal+box_i->nremote;
  }


  //define mapping between main force array and truncated version of threadprivate
  schedule->pfmap[tid] = ddcRealloc(schedule->pfmap[tid], schedule->lastdep[tid]*sizeof(int));
  int count = 0;
  for (int i = 0;i < schedule->lastdep[tid];i++)
  {
    OPT_BOX *box_tmp = boxlist + schedule->deplist[tid][i];
    schedule->pfmap[tid][i] = count;
    count+=box_tmp->nlocal+box_tmp->nremote;
  }
}

SCHEDULE* buildSchedule_MCAS(THREAD_SCHEDULER* this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize)
{
   SCHEDULE *schedule = schedule_init(boxSize, this->nThreads);
   int nbox = boxSize.nbox;
   int nThreads = this->nThreads;
   //double *boxtime = this->boxtime;
   CACHELINE_DOUBLE *boxtime = this->boxtime;

   int actual_nThreads = nThreads;


   int cnt_box = 0;
   for (int i=0; i<nbox; i++)
      if (boxlist[i].edge == 0)
	 cnt_box++;
	
   // To avoid infinite loop in first box assignment, nThreads must <= cnt_box.
   if (nThreads > cnt_box)
   {
      actual_nThreads = cnt_box;
   }
	
   schedule->nIdleThreads = nThreads - actual_nThreads;

   double workload[nThreads];
   for (int i=0;i<nThreads;i++)
   {
      workload[i] = 0.0;
      schedule->lastdep[i] = 0;
      schedule->lastalloc[i] = 0;
      schedule->pcount[i] = 0;
   }
	
   THREE_VECTOR ct[actual_nThreads];
   int clustersize[actual_nThreads];
   int box_owner[nbox];
   int changeOwnerStatus[nbox];

   int **growthlist = ddcMalloc(actual_nThreads*sizeof(int*));
   int *lastgrowth = ddcMalloc(actual_nThreads*sizeof(int));

   for (int i = 0;i < actual_nThreads;i++) {
      growthlist[i] = ddcMalloc(nbox*sizeof(int)); 
      lastgrowth[i] = 0;
   }

   for (int i = 0;i < nbox;i++)
   {
      box_owner[i]= -1;
      changeOwnerStatus[i] = 0;
   }
		
   // this loop selects an initial centroid for each thread.
   for (int k=0; k<actual_nThreads; k++)
   {  
      OPT_BOX *box;
      int rnd;
      do
      {
	 rnd=(rand() % nbox);
	 box = boxlist + rnd;
      } while (boxlist[rnd].edge == 1 || box_owner[rnd] > -1);
      //if (k==0) rnd = 1 + nx*(1 + ny);		
      box = boxlist+rnd;
      ct[k].x = box->x;
      ct[k].y = box->y;
      ct[k].z = box->z;
      box_owner[rnd] = k;
      workload[k] += boxtime[rnd].val; 	
      schedule->workalloc[k][schedule->lastalloc[k]++] = rnd;
      clustersize[k]=1;
      //get deplist

      ts_addWork(boxlist, growthlist, lastgrowth, rnd, k, boxSize, this->fullNN,this->profileThread); 
   }

/*
   if (getRank(0)==0)
   {
      printf("Seed boxes are:");   
      for (int i = 0;i < actual_nThreads;i++) 
         printf("%5d ",schedule->workalloc[i][0]);
      printf("\n");
   }
*/ 
	
   for (int mm = actual_nThreads;mm < cnt_box;mm++)
   {
      int minwork = 0;
      for (int k = 1;k < actual_nThreads;k++)
	 if (workload[k] < workload[minwork])
	    minwork = k;

      OPT_BOX *box; 
      double mindist= -1.0;	
      int minbox=0;
      
      for (int i=0;i < lastgrowth[minwork];i++)
      {
         if (box_owner[growthlist[minwork][i]] > -1 || boxlist[growthlist[minwork][i]].edge > 0) continue;	 
      	 box = boxlist + growthlist[minwork][i];	
	 double dx = ct[minwork].x-box->x;
	 double dy = ct[minwork].y-box->y;
	 double dz = ct[minwork].z-box->z;
	 double dist = dx*dx+dy*dy+dz*dz;
	 if (mindist > dist || mindist == -1.0)
	 {
	    minbox = growthlist[minwork][i];
	    mindist = dist;
	 }	
      }	
      if (mindist == -1.0)
      { //no available box to growth; randomly grab new one;
         if (this->strictSchedule==0)
         {
	    for (int i =0;i < nbox;i++)
	    {
	       if (box_owner[i]== -1 && boxlist[i].edge==0)
	       {
	          minbox = i;
	          i = nbox;
	       }		
	    }
	    //relocate the centroid to the new cluster
            box = boxlist + minbox;
	    ct[minwork].x = box->x;
	    ct[minwork].y = box->y;
	    ct[minwork].z = box->z;
	    clustersize[minwork]=0;
         }
         else //strictSchedule: no relocate centroid
         {
 	    mindist = -1.0;
            for (int i=0;i < lastgrowth[minwork];i++)
            {
	      if ((changeOwnerStatus[growthlist[minwork][i]] > 0) || (box_owner[growthlist[minwork][i]] == minwork) || (boxlist[growthlist[minwork][i]].edge > 0)) continue;	 
	      //if (box_owner[growthlist[minwork][i]] == minwork || boxlist[growthlist[minwork][i]].edge > 0) continue;	 
   	      box = boxlist + growthlist[minwork][i];	
	      double dx = ct[minwork].x-box->x;
	      double dy = ct[minwork].y-box->y;
	      double dz = ct[minwork].z-box->z;
	      double dist = dx*dx+dy*dy+dz*dz;
	      if (mindist > dist || mindist == -1.0)
	      {
	         minbox = growthlist[minwork][i];
	         mindist = dist;
	      }
            }
	      assert(mindist > -1.0);	
              changeOwnerStatus[minbox] = 1; 
	      int oldOwner = box_owner[minbox]; 
              assert(oldOwner!=minwork); 

		//	      printf("oldOwner oldworkload=%10.5le targetBox time = %10.5le\n",workload[oldOwner],boxtime[minbox]);
		
	      clustersize[oldOwner]--;	
      	      box = boxlist + minbox;
              //recalculate centroid for old box's owner thread
              ct[oldOwner].x = (ct[oldOwner].x*schedule->lastalloc[oldOwner]-box->x)/(clustersize[oldOwner]);  
              ct[oldOwner].y = (ct[oldOwner].y*schedule->lastalloc[oldOwner]-box->y)/(clustersize[oldOwner]);  
              ct[oldOwner].z = (ct[oldOwner].z*schedule->lastalloc[oldOwner]-box->z)/(clustersize[oldOwner]);  
	      
	      ts_removeWork(schedule,boxlist,growthlist,lastgrowth,&workload[oldOwner],boxtime[minbox].val,minbox,oldOwner,boxSize,this->fullNN); 
              
      	      mm--; //no progress for this iteration
	      //printf("at iteration %5d, box[%5d/%5d] change owner from tid[%5d] to tid[%5d]. Old owner workload = %10.5le\n",mm,minbox,nbox,oldOwner,minwork,workload[oldOwner]);
	 }	
      }	
   
      
   
		
      box = boxlist + minbox;
      //recalculate centroid
      ct[minwork].x = (ct[minwork].x*schedule->lastalloc[minwork]+box->x)/(1+clustersize[minwork]);  
      ct[minwork].y = (ct[minwork].y*schedule->lastalloc[minwork]+box->y)/(1+clustersize[minwork]);  
      ct[minwork].z = (ct[minwork].z*schedule->lastalloc[minwork]+box->z)/(1+clustersize[minwork]);  
		
      schedule->workalloc[minwork][schedule->lastalloc[minwork]++]=minbox;
      workload[minwork]+=boxtime[minbox].val;
		
      clustersize[minwork]++;
      box_owner[minbox]=minwork;
		
      ts_addWork(boxlist, growthlist, lastgrowth, minbox, minwork, boxSize, this->fullNN, this->profileThread); 
   } //end schedule

   //build deplist 
   ts_buildDeptList(schedule,boxlist,actual_nThreads,nbox);	
	
   //sort workalloc
   for (int k=0;k < actual_nThreads;k++)
   {
      qsort(schedule->workalloc[k], schedule->lastalloc[k], sizeof(int), intSortFunction);
   }

   for (int k=0;k < actual_nThreads;k++)
   	ddcFree(growthlist[k]);

   ddcFree(growthlist);
   ddcFree(lastgrowth);

   for (int k=0; k<nThreads; k++)
      schedule->workalloc[k][schedule->lastalloc[k]]= -1;
/*
   printf("Thread workload:");
   for (int k=0;k < actual_nThreads;k++)
   {  
	printf(" %10.5f",workload[k]);
   }
   printf("\n");

   printf("Each thread assignment:");
   for (int k=0;k < actual_nThreads;k++)
   {  
	printf(" tid = %d/%d lastalloc=%d",k,actual_nThreads,schedule->lastalloc[k]);

   	for (int kk=0;kk < schedule->lastalloc[k];kk++)
	{
		if (kk % 10 == 0) printf("\n\t");
		printf("%5d ",schedule->workalloc[k][kk]);
	}
	printf("\n");
   }
   printf("\n");
*/


   return schedule;
}

//Scheduling by Breath-First Allocation Scheduling (BFAS)
//The algorithm begin by randomly assign seed box for each thread. Put all neighbor boxes of the seed to queue of each thread.   
//Then, the least-load thread choose (and always choose) the head element in the queue.
//If the head element in the queue is already been picked by another thread, then choose the next one.
//Everytime the new box is assigned, put their neighbor boxes into a queue.
//If no element left in the queue, then pickup new unassigned box as a new seed.
//loop until all boxes are assigned.  
SCHEDULE* buildSchedule_BFAS(THREAD_SCHEDULER *this, OPT_BOX* boxlist, OPT_BOX_SIZE boxSize)
{
   SCHEDULE *schedule = schedule_init(boxSize, this->nThreads);
   int nbox = boxSize.nbox;
   int nThreads = this->nThreads;
   //double *boxtime = this->boxtime;
   CACHELINE_DOUBLE *boxtime = this->boxtime;

   int actual_nThreads = nThreads;

   int cnt_box = 0;
   for (int i=0; i<nbox; i++)
      if (boxlist[i].edge == 0)
	 cnt_box++;
	
   // To avoid infinite loop in first box assignment, nThreads must <= cnt_box.
   if (nThreads > cnt_box)
   {
      actual_nThreads = cnt_box;
   }
   schedule->nIdleThreads = nThreads - actual_nThreads;

   double workload[nThreads];
   for (int i=0;i<nThreads;i++)
   {
      workload[i] = 0.0;
      schedule->lastdep[i] = 0;
      schedule->lastalloc[i] = 0;
      schedule->pcount[i] = 0;
   }

	
   //this is counting non-edge box
   int queue_pnt[actual_nThreads];
   int box_owner[nbox];
   int changeOwnerStatus[nbox];
	
   int **growthlist = ddcMalloc(actual_nThreads*sizeof(int*));
   int *lastgrowth = ddcMalloc(actual_nThreads*sizeof(int));

   for (int i = 0;i < actual_nThreads;i++) {
      growthlist[i] = ddcMalloc(nbox*sizeof(int)); 
      lastgrowth[i] = 0;
   }

   for (int i = 0;i < nbox;i++)
   {
      box_owner[i]= -1;
      changeOwnerStatus[i] = 0;
   }

   // this loop selects an initial seed box for each thread.
   for (int k=0; k<actual_nThreads; k++)
   {  
      int rnd;
      do
      {
	 rnd=(rand() % nbox);
	 // OPT_BOX *box = boxlist + rnd;
      } while (boxlist[rnd].edge == 1 || box_owner[rnd] > -1);
      //if (k==0) rnd = 1 + nx*(1 + ny);		
		
      box_owner[rnd] = k;
      workload[k] += boxtime[rnd].val; 	
      schedule->workalloc[k][schedule->lastalloc[k]++] = rnd;
      //get deplist

      ts_addWork(boxlist, growthlist, lastgrowth, rnd, k, boxSize, this->fullNN, this->profileThread); 

      queue_pnt[k] = 0;
   }
/*	
   if (getRank(0)==0)
   {
      printf("Seed boxes are:");   
      for (int i = 0;i < actual_nThreads;i++) 
         printf("%5d ",schedule->workalloc[i][0]);
      printf("\n");
   }
*/ 
	
   for (int mm = actual_nThreads;mm < cnt_box;mm++)
   {
      int minwork = 0; //minimum workload threadid, initially set to thread 0
      int targetbox=-1; //index of target box to be added to least-loaded thread  
      for (int k = 1;k < actual_nThreads;k++)
	 if (workload[k] < workload[minwork])
	    minwork = k;

      while (queue_pnt[minwork] < lastgrowth[minwork])
      {
	 if (box_owner[growthlist[minwork][queue_pnt[minwork]]] == -1 && boxlist[growthlist[minwork][queue_pnt[minwork]]].edge==0) 
	    break;
	 else 
	    queue_pnt[minwork]++;
      }


      //if neighbor boxes is available, save index of to-be-added box    	
      if  (queue_pnt[minwork] < lastgrowth[minwork]) 
      {
	 targetbox = growthlist[minwork][queue_pnt[minwork]++];
      } 
      else //in case that all boxes from deplist[minwork] are already assigned or occupied by other threads, then relocate root of the thread  
      {
         if (this->strictSchedule==0)
         {
	    for (int i =0;i < nbox;i++)
	    {
	       if (box_owner[i]== -1 && boxlist[i].edge==0)
	       {
	          targetbox = i;
	          i = nbox;
	          queue_pnt[minwork]=lastgrowth[minwork]-1;
	       }		
   	    }
         }
         else
         {
            targetbox = -1;
            for (int i = 0;i < lastgrowth[minwork];i++)
            {
	       if (box_owner[growthlist[minwork][i]] != minwork  && boxlist[growthlist[minwork][i]].edge==0 && changeOwnerStatus[growthlist[minwork][i]] == 0) 
	       {
                  targetbox = growthlist[minwork][i];
                  i = lastgrowth[minwork];
               } 
	    else 
	       queue_pnt[minwork]++;
            }

            assert(targetbox >= 0);
            changeOwnerStatus[targetbox] = 1;
            int oldOwner = box_owner[targetbox];
            assert(oldOwner != minwork);
	    	
	    ts_removeWork(schedule,boxlist,growthlist,lastgrowth,&workload[oldOwner],boxtime[targetbox].val,targetbox,oldOwner,boxSize,this->fullNN);

            //reset queue for oldOwner;
            queue_pnt[oldOwner] = 0;
             
            mm--;
         }
      }
		
      schedule->workalloc[minwork][schedule->lastalloc[minwork]++]=targetbox;
      workload[minwork]+=boxtime[targetbox].val;
		
      box_owner[targetbox]=minwork;
		
      ts_addWork(boxlist, growthlist, lastgrowth, targetbox, minwork, boxSize, this->fullNN, this->profileThread); 
      } //end scheduling
	
   //build deplist 
   ts_buildDeptList(schedule,boxlist,actual_nThreads,nbox);	

   //sort workalloc
   for (int k=0;k < actual_nThreads;k++)
   {
      qsort(schedule->workalloc[k], schedule->lastalloc[k], sizeof(int), intSortFunction);
   }
	
   for (int k=0;k < actual_nThreads;k++)
   	ddcFree(growthlist[k]);

   ddcFree(growthlist);
   ddcFree(lastgrowth);


   for (int k=0; k<nThreads; k++)
      schedule->workalloc[k][schedule->lastalloc[k]]= -1;


   return schedule;
}


SCHEDULE* buildSchedule_trivial(THREAD_SCHEDULER* this,
				OPT_BOX* boxlist, OPT_BOX_SIZE boxSize)
{
   SCHEDULE* schedule = schedule_init(boxSize, this->nThreads);

   int nx = boxSize.nx;
   int ny = boxSize.ny;
   int nz = boxSize.nz;

   schedule->pcount[0] = 0;
   schedule->lastalloc[0] = 0;
        
   for (int iz=1; iz<nz-1; iz++)
   {
      int i =  1 + nx*(1 + ny *iz);
      if (boxlist[i].edge != 0 )
	 continue;
      for (int iy=1; iy<ny-1; iy++)
      {
	 for (int ix=1; ix<nx-1; ix++)
	 {
	    OPT_BOX *box_i = boxlist+i;
	    int nTotal_i = box_i->nlocal+box_i->nremote;
	    if (nTotal_i <= 0)
	    {
	       i++;
	       continue;
	    }
	    schedule->workalloc[0][schedule->lastalloc[0]++]=i;
	    i++;
	 }
	 i+=2;
      }
   }

   schedule->workalloc[0][schedule->lastalloc[0]] = -1;
        
   for (int i=0; i<boxSize.nbox; i++)
      schedule->deplist[0][i] = i;      
   schedule->lastdep[0] = boxSize.nbox;
   schedule->nIdleThreads = 0;	

   return schedule;
}

void ts_addWork(OPT_BOX *boxlist, int **growthlist, int *lastgrowth, int targetBoxId, int targetThread, OPT_BOX_SIZE boxSize, int fullNN, int profileThread) 
{
      int nntable[27];
      int nntablecount = 0;
      OPT_BOX *box = boxlist+targetBoxId;	
      nntable[nntablecount++]=targetBoxId;

      if (fullNN == 1) 
      {	
         int ix = box->x;	
         int iy = box->y;	
         int iz = box->z;	

	 int nx = boxSize.nx;
	 int ny = boxSize.ny;
	 int nz = boxSize.nz;
         
  	   for(int z= -1;z<=1;z++) 
	   for(int y= -1;y<=1;y++) 
	   for(int x= -1;x<=1;x++) 
	   {
		int l = x + nx*(y + ny*z);
		if (l == 0 )  continue; 
		int jx = ix + x; 
		int jy = iy + y; 
		int jz = iz + z; 
		if (jx < 0 || jx > nx-1 || jy < 0 || jy > ny -1 || jz < 0 || jz > nz-1 ) continue ; 
		int j = jx + nx*(jy + ny*jz);
	 	nntable[nntablecount++]=j;
           }
      }
      else
      {

         for  (int j=0;j <box->nn;j++ ) /* loop over neighboring boxes of box */
         {
	    nntable[nntablecount++]=targetBoxId+box->nlist[j].offset;
         }
      }
		
        // add to growthlist
      ts_addGrowthList(growthlist[targetThread],nntable,&lastgrowth[targetThread],nntablecount);	

      if (profileThread == 1)
      {
         fprintf(profileFile,"%d %d %d %d\n",targetThread,box->x,box->y,box->z);
      }	
}

void ts_removeWork(SCHEDULE *sc, OPT_BOX *boxlist, int **growthlist, int *lastgrowth, double *threadWorkload, double targetBoxtime, int targetBoxId, int targetThread, OPT_BOX_SIZE boxSize, int fullNN) 
{
   *threadWorkload -= targetBoxtime;
   int isFound = 0;
   for (int i = 0;i < sc->lastalloc[targetThread];i++)
   {
      if (sc->workalloc[targetThread][i] == targetBoxId)    
      {
         sc->lastalloc[targetThread]--;
         sc->workalloc[targetThread][i] = sc->workalloc[targetThread][sc->lastalloc[targetThread]];
         i = sc->lastalloc[targetThread];
	 isFound = 1;
      }
   }     	
   assert(isFound == 1);
   // rebuild thread growth list
   ts_buildThreadGrowthList(sc, boxlist, growthlist[targetThread], &lastgrowth[targetThread], targetThread, boxSize, fullNN);

}

void ts_addGrowthList(int *growthlist,int *nntable,int *lastgrowth, int nntablecount)
{
      int lastgrowth_s = *lastgrowth;	
      for (int i = 0;i < nntablecount;i++)
      {
	 int isfound = 0;
	 for (int j = 0;j < lastgrowth_s;j++)
	 {
	    if (nntable[i]==growthlist[j])
	    {
	       isfound = 1;
	       j = lastgrowth_s;
	    }
	 }
			
	 if (isfound == 0)
	    growthlist[lastgrowth_s++]=nntable[i];
      }
      *lastgrowth = lastgrowth_s;	 
}

void ts_buildThreadGrowthList(SCHEDULE *sc, OPT_BOX *boxlist, int *growthlist, int *lastgrowth, int targetThreadId,OPT_BOX_SIZE boxSize, int fullNN)
{  
   /*
   depStatus = 0 -> not in growthlist and not in workalloc 
   depStatus = 1 ->     in growthlist but not in workalloc 
   depStatus = 2 ->     in growthlist and     in workalloc 
   */

   int deptStatus[boxSize.nbox];

   int *workallocThread = sc->workalloc[targetThreadId]; 


   for (int i=0;i < boxSize.nbox;i++)
      deptStatus[i] = 0;

   int lastalloc_s = sc->lastalloc[targetThreadId];
   for (int i = 0;i < lastalloc_s;i++)
   {
      int kk = workallocThread[i];
      OPT_BOX *box = boxlist+kk;
      //deptStatus[kk] = 2; //
      deptStatus[kk] = 2; //
      if (fullNN == 1) 
      {	
         int ix = box->x;	
         int iy = box->y;	
         int iz = box->z;	

         int nx = boxSize.nx;
         int ny = boxSize.ny;
	 int nz = boxSize.nz;
         
  	 for(int z= -1;z<=1;z++) 
	 for(int y= -1;y<=1;y++) 
	 for(int x= -1;x<=1;x++) 
	 {
	    int l = x + nx*(y + ny*z);
	    if (l == 0 )  continue; 
	    int jx = ix + x; 
	    int jy = iy + y; 
	    int jz = iz + z; 
	    if (jx < 0 || jx > nx-1 || jy < 0 || jy > ny -1 || jz < 0 || jz > nz-1 ) continue ; 
	    int j = jx + nx*(jy + ny*jz);
	    if (deptStatus[j] == 0) deptStatus[j] = 1;
         }
      }
      else
      {
         for  (int j=0;j <box->nn;j++ ) /* loop over neighboring boxes of box */
         {
            int ix = kk+box->nlist[j].offset;
	    if (deptStatus[ix] == 0) deptStatus[ix] = 1;
         }
      }


      int lastgrowth_s = 0;
      for (int i = 0;i < boxSize.nbox;i++)
      {
         if (deptStatus[i] == 1) growthlist[lastgrowth_s++] = i;
      }
      *lastgrowth = lastgrowth_s;
   }   
}

void ts_buildDeptList(SCHEDULE *sc, OPT_BOX *boxlist, int actual_nThreads, int nbox)	
{
   int *deptStatus = ddcMalloc(nbox*sizeof(int));

   for (int k=0;k < actual_nThreads;k++)
   {
      int *workallocThread = sc->workalloc[k]; 

      for (int i=0;i < nbox;i++)
         deptStatus[i] = 0;

      int lastalloc_s = sc->lastalloc[k];
      for (int i = 0;i < lastalloc_s;i++)
      {
	 int ix = workallocThread[i];
         OPT_BOX *box = boxlist+ix;
         deptStatus[ix] = 1;
         for  (int j=0;j <box->nn;j++ ) /* loop over neighboring boxes of box */
         {
            deptStatus[ix+box->nlist[j].offset] = 1;
         }
      } 

      for (int i = 0;i < nbox;i++)
      {
         if (deptStatus[i] > 0) sc->deplist[k][sc->lastdep[k]++] = i;
      }
   }

   ddcFree(deptStatus);
}


