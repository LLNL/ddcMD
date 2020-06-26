#ifndef PRIMES_H
#define PRIMES_H

/** This is a prime number generator.  Because all state is stored
 *  internally in static variables, it is not thread safe.  However, it
 *  can be called by any combination of clients that need prime numbers.
 *  There is no MPI communication in the generator or the
 *  initialization.
 *
 *  These routines are written so that each task on the system has a set
 *  of ranges of integers (or blocks) that will only be searched by that
 *  task.  If the task exhausts all of the primes in a block it moves on
 *  to the next block assigned to that task.  All of this occurs without
 *  any communication between the tasks.  Hence it is safe to call
 *  nextPrime() on any task at any time (after prime_init).
 *  
 *  prime_init needs to be called by all tasks that will need prime
 *  numbers.  taskId needs to be a unique number in the range
 *  0<=taskId<nTasks and nTasks is the number of tasks that will
 *  generate primes.  It is safest to call prime_init on all tasks even
 *  if there are some that won't generate primes.
 */
void prime_init(unsigned blockSize, unsigned taskId, unsigned nTasks);
unsigned long long nextPrime(void); 

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
