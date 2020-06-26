/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <sys/time.h>

#include "rectimer.h"

#include "genptr.h"

#ifdef __bgq__
  #include <hwi/include/bqc/A2_inlines.h>
  typedef long long int timerval_t;

  timerval_t gettime(void) {
    return GetTimeBase();
  }

  timerval_t gettime_tick(void) {
    const timerval_t million = 1000000,f = 1600;
    return f*million;
  }
#elif 1
  /* #error "Supposed to hit __bgq__ define , not here" */
  #include <sys/time.h>
  typedef double timerval_t;

  double gettime(void) {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
  }

  double gettime_tick(void) {
    return 1.0;
  }
#else
  #include <time.h>
  typedef long long int timerval_t;

  timerval_t gettime(void) {
    return clock();
  }

  timerval_t gettime_tick(void) {
    return CLOCKS_PER_SEC;
  }
#endif

typedef struct {
  int itemsize,size,nalloc,used;
  void **key;
  int *next;
  void *data;
  int *hashtable;
} voidhash;

typedef struct {
  double ttot;
  int nstart,nstop;
} callnode;

typedef struct REC_TIMER_S {
  const char *name;
  int nstart,nstop;
  timerval_t ttot;
  voidhash callgraph;
} timer;

typedef struct {
  thandle th;
  timerval_t starttime;
} stackitem;

static struct {
  int nused,nalloc;
  thandle *list;
} timers = { 0 , 0 , NULL };

static struct {
  int nused,nalloc;
  stackitem *list;
} timerstack = { 0 , 0 , NULL };

static int isprime(int x) {
  int k = 3;
  if((x&1) == 0) return 0;
  while(k*k <= x) {
    if(x % k == 0) return 0;
    k += 2;
  }
  return 1;
}
static int nextprime(int x) {
  while(!isprime(x)) x++;
  return x;
}
static void voidhash_new(voidhash *vh,int itemsize,int nitem) {
  int i;
  vh->itemsize = itemsize;
  vh->size = nextprime(nitem);
  vh->nalloc = 0;
  vh->used = 0;
  vh->key = NULL;
  vh->next = NULL;
  vh->data = NULL;
  vh->hashtable = (int *) malloc(sizeof(int) * vh->size);  
  for(i = 0; i<vh->size; i++)
    vh->hashtable[i] = -1;
}
/*
static void void_hash_delete(voidhash * vh) {
  free(vh->hashtable);
  free(vh->key);
  free(vh->next);
  free(vh->data);
  vh->itemsize = vh->size = vh->nalloc = vh->used = 0;
  vh->hashtable = NULL;
  vh->key = NULL;
  vh->next = NULL;
  vh->data = NULL;
}
*/
static void * voidhash_lookup_insert(voidhash *vh,void *key,void *initializer) {
  int idx = ((size_t) key) % vh->size;
  int p = vh->hashtable[idx];
  while(p >= 0 && vh->key[p] != key)
    p = vh->next[p];
  if(p < 0) {
    if(vh->used >= vh->nalloc) {
      int n = vh->nalloc * 2;
      if(n < 10) n = 10;
      vh->key = (void **) realloc(vh->key,sizeof(void *) * n);
      vh->next = (int *) realloc(vh->next,sizeof(int) * n);
      vh->data = (void *) realloc(vh->data,vh->itemsize * n);
      vh->nalloc = n;
    }
    p = vh->used++;
    vh->key[p] = key;
    vh->next[p] = vh->hashtable[idx];
    vh->hashtable[idx] = p;
    memcpy(idx_ptr(vh->data,p,vh->itemsize),initializer,vh->itemsize);
  }
  return idx_ptr(vh->data,p,vh->itemsize);
}


thandle timer_lookup_insert(const char *name) {
  int i;
  thandle th;
  for(i = 0; i<timers.nused; i++)
    if(strcmp(timers.list[i]->name,name) == 0) break;

  if(i < timers.nused)
    th = timers.list[i];
  else {
    if(timers.nused >= timers.nalloc) {
      int n = timers.nalloc * 2;
      if(n < 10) n = 10;
      timers.list = (thandle *) realloc(timers.list,sizeof(*timers.list) * n);
      timers.nalloc = n;
    }
    th = (thandle) malloc(sizeof(timer));
    timers.list[timers.nused++] = th;
    th->nstart = th->nstop = 0;
    if(0) {
      th->name = strdup(name);
    } else {
      char *s = (char *) malloc((strlen(name)+1)*sizeof(char));
      strcpy(s,name);
      th->name = s;
    }
    th->ttot = 0;
    voidhash_new(&(th->callgraph),sizeof(callnode),13);
  }
  return th;
}

thandle rectimer_start(thandle th,const char *name) {
  if(th == NULL) th = timer_lookup_insert(name);
  th->nstart++;
  if(timerstack.nused >= timerstack.nalloc) {
    timerstack.nalloc *= 2;
    if(timerstack.nalloc < 10) timerstack.nalloc = 10;
    timerstack.list = (stackitem *) realloc(timerstack.list,sizeof(*timerstack.list) * timerstack.nalloc);
  }

  if(timerstack.nused > 0) {
    /* Accumulate child time in call graph of parent timer */
    callnode zero = { 0.0 , 0 , 0 };
    callnode *node = (callnode *)
      voidhash_lookup_insert(&(timerstack.list[timerstack.nused-1].th->callgraph),th,&zero);
    node->nstart++;
  }

  timerstack.list[timerstack.nused].th = th;
  timerstack.list[timerstack.nused++].starttime = gettime();

  return th;
}

void rectimer_stop(thandle th) {
  const int depth = timerstack.nused - 1;
  timerval_t tnow = gettime();
  timerval_t tcall = tnow - timerstack.list[depth].starttime;
  assert(th != NULL);
  assert(th == timerstack.list[depth].th);
  th->nstop++;
  if(th->nstart == th->nstop)
    th->ttot += tcall;
  timerstack.nused = depth;

  if(depth > 0) {
    /* Accumulate child time in call graph of parent timer */
    /*
    callnode zero = {0.0 , 0};
    callnode *node = (callnode *)
      voidhash_lookup_insert(&(timerstack.list[depth-1].th->callgraph),th,&zero);
    */
    /*
    printf("stop: %s.%s updated from (%.3e,%d) ",
	   timerstack.list[depth-1].th->name,th->name,node->ttot,node->ncalls);
    */
    {
      callnode *node = (callnode *)
	voidhash_lookup_insert(&(timerstack.list[depth-1].th->callgraph),th,NULL);
      node->nstop++;
      if(node->nstart == node->nstop)
	node->ttot += tcall;
    }
    //printf("to (%.3e,%d)\n",node->ttot,node->ncalls);
  }
}

int rectimer_ncalls(thandle th) { assert(th != NULL); return th->nstop; }
double rectimer_ttot(thandle th) {
  assert(th != NULL);
  return th->ttot / (double) gettime_tick();
}
double rectimer_tavg(thandle th) {
  assert(th != NULL);
  if(th->nstop > 0) {
    return ((double) th->ttot)/th->nstop / (double) gettime_tick();
  }
  return 0.0;
}

void rectimer_printreport(FILE *fp) {
  int i;
  fprintf(fp,"--- TIMING REPORT ---\n");
  for(i = 0; i<timers.nused; i++) {
    int j;
    thandle th = timers.list[i];
    fprintf(fp,
	    "  %-24s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
	    th->name,rectimer_ttot(th),rectimer_tavg(th),(double) rectimer_ncalls(th));
    for(j = 0; j<th->callgraph.used; j++) {
      thandle ch = (thandle) th->callgraph.key[j];
      callnode *node =  (callnode *) idx_ptr(th->callgraph.data,j,sizeof(callnode));
      int nj = node->nstop;
      double tot = node->ttot / (double) gettime_tick();
      double avg = 0.0;
      if(nj > 0) avg = tot/nj;
      fprintf(fp,
	      "    %-22s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
	      ch->name,tot,avg,(double) nj);
    }
    if(j > 0)
      fprintf(fp,"\n");
  }
}


char * rectimer_printreport_string(void) {
  int ptr = 0,len = 0,flag;
  char *str = NULL;
  int i;

  for(flag = 0; flag<=1; flag++) {
    ptr += snprintf(str+ptr,flag*(len-ptr),"--- TIMING REPORT ---\n");
  
    for(i = 0; i<timers.nused; i++) {
      int j;
      thandle th = timers.list[i];
      ptr += snprintf(str+ptr,flag*(len-ptr),
		      "  %-24s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
		      th->name,rectimer_ttot(th),rectimer_tavg(th),(double) rectimer_ncalls(th));
      for(j = 0; j<th->callgraph.used; j++) {
	thandle ch = (thandle) th->callgraph.key[j];
	callnode *node =  (callnode *) idx_ptr(th->callgraph.data,j,sizeof(callnode));
	int nj = node->nstop;
	double tot = node->ttot / (double) gettime_tick();
	double avg = 0.0;
	if(nj > 0) avg = tot/nj;
	ptr += snprintf(str+ptr,flag*(len-ptr),
			"    %-22s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
			ch->name,tot,avg,(double) nj);
      }
      if(j > 0)
	ptr += snprintf(str+ptr,flag*(len-ptr),"\n");
    }
    if(flag == 0) {
      len = ptr+1;
      str = (char *) malloc(len * sizeof(char));
      ptr = 0;
    }
  }
  return str;
}


void rectimer_printreport_file(fprintf_fun_type fprintf_fun,void *fp) {
  int i;
  fprintf_fun(fp,"--- TIMING REPORT ---\n");
  for(i = 0; i<timers.nused; i++) {
    int j;
    thandle th = timers.list[i];
    fprintf_fun(fp,
	    "  %-24s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
	    th->name,rectimer_ttot(th),rectimer_tavg(th),(double) rectimer_ncalls(th));
    for(j = 0; j<th->callgraph.used; j++) {
      thandle ch = (thandle) th->callgraph.key[j];
      callnode *node =  (callnode *) idx_ptr(th->callgraph.data,j,sizeof(callnode));
      int nj = node->nstop;
      double tot = node->ttot / (double) gettime_tick();
      double avg = 0.0;
      if(nj > 0) avg = tot/nj;
      fprintf_fun(fp,
	      "    %-22s :  tot = %15.5e  avg = %15.5e  ncall = %15.5e\n",
	      ch->name,tot,avg,(double) nj);
    }
    if(j > 0)
      fprintf_fun(fp,"\n");
  }
}
