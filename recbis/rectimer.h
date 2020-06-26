/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef RECTIMER__
#define RECTIMER__

#include <stdio.h>

typedef struct REC_TIMER_S *thandle;

thandle rectimer_start(thandle th,const char *name);
/*
   Start a timer with name 'name'. If th is NULL a timer handle is
   returned. If th is not NULL, th is retuned name is ignored. th
   is then assumed to be already initialized (i.e. a return value
   from a previous call to rectimer_start()).

   Recursive calls are supported, and recursion statistics for each
   timer is provided through the report printing  functions below.

   The file 'functimer.h' defines the macros STARTTIMER and STOPTIMER
   that can be put at the top and before returning in a function to
   record total time spent in that function with minial code addition.
*/
void rectimer_stop(thandle th);
/* Stop the timer designated by th. */

int rectimer_ncalls(thandle th);
/* Return how many times th has been stopped */

double rectimer_ttot(thandle th);
/* Return total amount of wall clock time th has been active (between start and stop) */

double rectimer_tavg(thandle th);
/* Return average amount of time th was active in a start-stop interval */

void rectimer_printreport(FILE *fp);
/* Print a report of timer statistics to a text file */

char * rectimer_printreport_string(void);
/* Return a string which contains a report of timer statistics */

typedef int (*fprintf_fun_type)(void *fileobject,const void *fmt,...);
void rectimer_printreport_file(fprintf_fun_type fprintf_fun,void *fp);
/*
   Print a report of timer statiscts to object fp, assuming fprintf_fun is
   a variable argument function that is called like fprintf. fp will be
   the first argument in each call to fprintf_fun. This allows printing
   reports to other types of file-like objects
*/

#endif
