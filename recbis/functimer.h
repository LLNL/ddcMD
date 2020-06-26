/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef FUNCTIMER__
#define FUNCTIMER__

/*
   Convenient macros for timing whole functions. Put STARTTIMER
   at the top of a function, and STOPTIMER just before each
   return statement. If you need additional timers in a function
   you can use the interface in rectimer.h to define any number
   of named timers.

   The report printing functions declared in rectimer.h can
   print statiscts on wall clock time that each timer has
   been active. Nesting and recursion is allowed, and discernable
   from the printed reports.
*/

#define STARTTIMER \
  static thandle t_all = NULL; \
  t_all = rectimer_start(t_all,__func__);
#define STOPTIMER rectimer_stop(t_all);

//#define STARTTIMER
//#define STOPTIMER
#endif
