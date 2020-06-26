#ifndef RANDOM_H
#define RANDOM_H
enum RANDOM_CLASS { URANDOM, LCG64 };
enum RANDOM_STATETYPE {SCALAR, REGISTERED};
#include <stdio.h>
#include "gid.h"
#include "object.h"
#include "three_algebra.h"

/** This is a long comment, but it contains important information.  If
 *  you're writing code that uses the RANDOM class, do yourself a favor
 *  and take the time to read it.
 *
 *
 *  RANDOM is a bit of a tricky class because it has to support a
 *  somewhat complex usage model.  RANDOM is structured as a
 *  "polymorphic" structure so that we could, in principle, support
 *  various random number generator implementation all through a single
 *  interface.  Usually in a polymorphic structure we bury all of the
 *  subclass parameters in a void* data member, typically named parms.
 *  The parms are passed to the subclass with each function call since
 *  the polymorphic functions all take a pointer to a base class as an
 *  argument (i.e., f(RANDOM* random) ) Within the subclass the void*
 *  gets cast to a subclass parms type and the subclass can pull out
 *  whatever data is needed.  In this method the caller of the base
 *  class never even needs to know that the void* parms exist.
 *
 *  RANDOM doesn't work that way because of the need to have an
 *  independent random number generator for each particle in the system.
 *  It might seem like overkill to have a separate random generator for
 *  each particle, however one of the key advantages is that it helps
 *  maintain reproducibility even when the code is run on different
 *  numbers of tasks.  Because the random generators are tethered to the
 *  particles, a particle always experiences the same random stream
 *  regardless of which task it is on or how many other particles are on
 *  the task.  This feature alone makes it worth the trouble.
 *
 *  Since we want to have one generator per particle, the first though
 *  would be to create an array of RANDOMs.  However, this isn't a great
 *  idea for two reasons.  First, with the overhead of the OBJECT stuff
 *  at the top of the struct, and all of the "virtual" functions,
 *  sizeof(RANDOM) isn't a small number so one RANDOM per particle would
 *  chew up quite a bit memory for large particle counts.  Furthermore,
 *  since we are going to move the RANDOM with the particle when it gets
 *  reassigned to another we want to minimize the amount of data we need
 *  to move.  (And we haven't even mentioned the complications that come
 *  from moving a variable size class with data that isn't contiguous in
 *  memory, including function pointers that only work if the address
 *  space is the same on all tasks.  But I digress.)
 *
 *  The solution to this problem is to store on a per particle basis,
 *  only the minimum amount of data that is needed.  That is the parms
 *  of the subclass.  We also modify our polymorphic function signatures
 *  so that they take a parms instead of (or in some cases in addition
 *  to) RANDOM*.  This explains why, for example, uniformRandom takes
 *  void* parms as it argument instead of the RANDOM* you might expect. 
 *
 *  So far, so good, but this introduces new complications.  First, what
 *  class will own the array of parms, one parm per particle?  It is
 *  tempting to put it in STATE, but we reject that approach because it
 *  isn't extensible.  What if we need two different generators per
 *  particle?  Instead we store the array of parms in the parmsArray
 *  member of the random class, and also take care of calling
 *  particleRegisterinfo (when random_registerParms is called) so that
 *  the parms will move with the particles.  One nice thing about this
 *  approach is that the random parms only get registered if some caller
 *  is actually using them.
 *
 *  IMPORTANT NOTE: calling random_registerParms does not initialize the
 *  parms.  The caller must initialize the parms through some other
 *  method such as calling parse, bread, defaultValue for each
 *  particle's parms.
 *
 *  ANOTHER IMPORTANT NOTE: While calling random_registerParms does
 *  result in a calls to particleRegisterinfo, is does not allocate
 *  memory.  There will be no memory available to store and retrieve
 *  data in the parmsArray until particleAllocateinfo is called.
 *  (particleAllocateinfo is typically invoked by calling resize.)
 *
 *  Second, how do callers get access to the parms for the ith particle?
 *  There is no point in trying to export the parmsArray to callers
 *  because it is a void*.  Pointer arithmatic on a void* is undefined,
 *  so it is impossible to know where the data for the ith particle is
 *  without knowing something about the concrete subclass type.  This
 *  would defeat the whole point of polymorphism.  Instead, we provide
 *  the function void* random_getParms(RANDOM*, unsigned i) to get the
 *  ith parms.  The pointer returned by this function can then be passed
 *  to uniformRandom, parse, bwrite, defaultValue, etc.
 *
 *  Third, what if a caller really just wants a single random number
 *  generator and not one per particle?  And don't forget that due to
 *  the way object.data is set up it is entirely likely that multiple
 *  objects all specified the same random object and some of them are
 *  particle based and some are not.  So, RANDOM must be prepared to do
 *  double duty as both a single stream and a particle matched
 *  generator, perhaps even simultaneously.  To support single stream
 *  usage we provide the function void* random_newParms(RANDOM*).
 *  Callers that need a single stream random should acquire the user
 *  requested RANDOM* from SYSTEM, then call random_newParms to get a
 *  single stream parms.  If there is some value for the parms that
 *  needs to be restored (as from a restart file) the caller can call
 *  parse (for example) to put the desired values in parms.  The single
 *  stream caller needs to keep track of both a RANDOM* random and a
 *  void* parms separately, (so that function calls like
 *  random->uniformRandom(parms) can be made) but this is the cost of
 *  doing business with a class that has to support a complex usage
 *  model.
 *
 *  Fourth, what if you want randomness that will give the same results
 *  regardless of the number of tasks the code is using, but you don't
 *  want to create a full blown set of particle based random objects.
 *  For example, you want to shuffle gids or set the velocities to a new
 *  Boltzmann distribution.  In this case you can use prand48.  The idea
 *  here is that each particle re-seeds the generator with a new seed
 *  that is based on its gid as well as a user supplied (or time of day
 *  generated) seed.  That way the particle always gets the same random
 *  values regardless of how many tasks are in operation.  For a good
 *  example of how this works see the thermalize transform.
 
 *  Finally, sometimes we just want a disposable generator to run off a
 *  few random numbers and then go away.  We often use drand48, but this
 *  is a bit of a hassle if we want different streams on each task, and
 *  we don't want to generate a RANDOM for each particle because the
 *  overhead is too high.  Furthermore, we may be in a context where the
 *  user hasn't provided the name of a RANDOM object for us to use.  In
 *  this case, if we don't care about repeatability on varied numbers of
 *  tasks we can call random_getDisposableRandomGenerator.  Note that
 *  there is no way to specify what type of generator you get.  We can
 *  add something to the argument list if that is ever important.  Once
 *  you're done, call random_free.
 *
 */

#ifdef __cplusplus
extern "C"{ 
#endif


typedef struct random_st
{
	char *name;
	char *objclass;
	char *value;
	char *type;		/* model */
	void *parent; 
	enum RANDOM_CLASS itype;	/* integer label for type */
	LONG64 seed;
   int useDefault; 
	unsigned sizeOfParms;
	void* parmsArray; 
	char *(*parse) (char *string, void *parms);
	void (*defaultValue) (unsigned n, LONG64 seed, void *parms);
	char *(*write) (void *parms);
	unsigned (*bread)(unsigned char* in, void *parm);
	unsigned (*bwrite)(unsigned char* out, void *parm);
//	void (*write_dynamics) (struct random_st *random, FILE *file);
	double (*uniformRandom)(void *parm); 
	TWO_VECTOR (*uniformRandom2)(void *parm); 
} RANDOM;

RANDOM *random_init(void *parent, char *name);
void random_registerParms(RANDOM* random);
void* random_newParms(RANDOM* random);

void *random_getParms(RANDOM *random, unsigned k);
int random_checkValue(RANDOM *random, unsigned k);

void random_getDisposableGenerator(RANDOM* self, void* parms);
void random_free(RANDOM* self, void* parms);

double       gasdev(RANDOM *random,void *parm);
THREE_VECTOR gasdev3d(RANDOM *random, void *parm);

THREE_VECTOR jutdev(RANDOM *random, void *parm, double Trel);

LONG64 generateRandomSeed(void);

typedef struct PRAND48_STATE_st
{
	unsigned short seedVec[3];
} PRAND48_STATE;

double       gasdev0(PRAND48_STATE*);
THREE_VECTOR gasdev3d0(double sigma, PRAND48_STATE*);

THREE_VECTOR jutdev0(PRAND48_STATE*, double Trel);

double prand48(PRAND48_STATE* self);
PRAND48_STATE prand48_init(LONG64 gid, LONG64 seed, LONG64 callSite);

void missingRandomError(char *string);
void testForRandomKeyword(OBJECT  *obj);

#ifdef __cplusplus
}
#endif


#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */

/* Here are some more lengthy notes on the refactoring of random. */
/* I'm putting them here so that I don't loose them.  Eventually we */
/* will get this all resolved and this stuff can be deleted or */
/* incorporated into the comments that explain the RANDOM class. */


/* Problem: Multiple groups all specify the same random object */
/* * random_init gets called by each group */
/*   - each group gets an independent instance of RANDOM* even though  */
/*     there is only one object in object.data */
/*   - This happens even when you have two different Langevin groups. */
/* * ownership of the per-particle random parms is unclear. */
/*   - if RANDOM owns the state then each independent instance of */
/*     RANDOM potentially adds a complete set of parms to the  */
/*     particle registry. */
/*   - if STATE owns the random parms then multiple groups could  */
/*     attempt to read/write the same data */
/*   - and what about instances of RANDOM that don't need */
/*     per-particle parms? */
/*   - and what if different groups want to use different types */
/*     of RANDOMs with different sized parms?  How do we put that */
/*     into STATE? */



/* Solution: */
/* * Add an array of RANDOM to the SYSTEM */
/* * Add random keyword to SYSTEM object.  This is a list of all the */
/*   RANDOMs that are in use. */
/*   - To avoid breaking every input deck that exists system_init will also */
/*     scan through all the groups it knows about and add the values of  */
/*     their random keywords to the list. */
/* * SYSTEM calls random_init for each RANDOM in the object.data. */
/* * GROUPS now need to get their RANDOMs from the SYSTEM instead of */
/*   by calling random_init. */
/*   - Question:  How do GROUPS ask for a random?  It must be something */
/*     like the way that COLLECTION asks for SPECIES or a GROUP. */




/* Remaining issues:  Who handles reading or writing of the random parms? */
/* should it be GROUP or RANDOM? */
/* * We don't need multiple groups to save the same data. */
/* * We don't want to initialize random parms for particles that never */
/*   use random numbers. */
/* * We don't want RANDOM to save parms for particles that aren't in use */
/*   (in variable records).  In fixed records we're going to pad. */
/* * If RANDOM does the reading how does it tell which states to "turn on" */

/* * Another thing: particleAllocateinfo gets memory from ddcMalloc.  So, */
/*   the RANDOM parms data will be uninitialized, hence there is no  */
/*   practical way to look at a given state and tell whether it has been */
/*   set up or not. */
/* * Initialization order matters.  Initialize RANDOM before GROUPS */

/* * It looks like Jim may have started to write some code to handle */
/*   initializing a scalar RANDOM from object.data.  random_init contains  */
/*   an object_get on "state".  However the code is clearly incomplete */
/*   (nothing is done with the string that is read) and there is no */
/*   corresponding write_dynamics code.  What to do with it?  Putting */
/*   a single state in the RANDOM object.data doesn't make sense unless */
/*   we're going to use a rule that a RANDOM can be either scalar or */
/*   particle but not both.  Then users will have to use a separately  */
/*   named RANDOM for the single stream RANDOMS.  Is that what we want? */

/* * Usually with polymorphic classes the subclass parms file is defined */
/*   only in the subclass .c file, not in the .h  It is after all, no */
/*   one else's business.  Should we do this in lcg64? */



/* Notes: */

/* 1.  Why isn't the default for system type NORMAL?  Why does system even */
/*     have a type?  What possible other values are imagined? */
/* 2.  system->changed appears to be unused except for being set once in  */
/*     system_init. */
/* 3.  Both group_init and species_init call object_initialize on the name */
/*     they are given.  This means the code will not fail if the user */
/*     fails to define a group or species that is named in SYSTEM. Instead */
/*     the code will silently create a default group or species.  This */
/*     isn't a very good way to check for errors. */
	
/* 6.  Do we always free the memory that is returned from object_getv? */
/*     Even when it isn't a string array? */

/* 7.  It looks like we might not need the nMember data in groups or */
/*     species any more.  Jim might be using them for something else */
/*     so check. */

/* 8.  Back in r1175 Jim created functions like group_bread, presumably to */
/*     be default values for members like group->bread.  However, it looks */
/*     like it was never actually hooked up.  Figure out which way to do */
/*     it.  The comments that I wrote awhile ago in group.h says that the  */
/*     functions are optionally defined by the subclasses, but it doesn't */
/*     say that the function pointers may be left NULL and need to be  */
/*     check before calling. */
    
/* 9.  radiation.c defines parse, write, bwrite, and bread functions even */
/*     though it apparently has no use for them.  I'm going to delete */
/*     them. */

/* 10. radiation.c appears to use the single stream random option. */
/*     It also tries to save the parms in the restart file using the  */
/*     write_dynamics method.  What can we do about the fact that unless  */
/*     we scheme otherwise every task will have a different parms?  We */
/*     clearly don't want to write all of them in the restart file.  Do  */
/*     we use the single stream on one task only?  What if we don't? */

/* 11. When radiation_parms intializes parms->randomParms from the */
/*     object database there is no error checking, including any check */
/*     that the random parms still match the type of random generator */
/*     that the group is using.  And what if the user tries to supply  */
/*     randomParms and messes it up? */

/* 12. lcg64_parse appears to set no error flag.  Isn't it supposed to?  */
/*     How else can we tell that the parms we are trying to parse is bogus? */

/* 13. There is some obviously wrong code in ionization_io.c involving a */
/*     call to lcg64_default.  It is preprocessed out, but I haven't take */
/*     time to figure out what is really going on. */

/* 14. Is langevin_default a correct implementation?  In particular, is */
/*     it a good idea to mix the seed with the label? */

/* 15. Why is there a langevin_write_dynamics function?  I can think of no */
/*     good purpose for it.  I deleted the unused random from the function. */
/*     I think we should delete the whole function. */

/* 16. Should we delete the first argument (unsigned n) from the */
/*     RANDOM::defaultValue function.  I think it no longer serves any */
/*     purpose. */

/* 17. There is no error check in system_getRandomByName.  Callers should */
/*     make sure that they get a non-NULL pointer back. */

/* 18. The signature of system_getRandomByName isn't very satisfactory. */
/*     The caller likely has to pull SYSTEM out of thin air.  Of course  */
/*     the alternatives (such as are used for GROUP and SPECIES) aren't */
/*     great either.  What to do.... */

/* 19. I got rid of the call to group->defaultValue(-1,...) in */
/*     collection_init.  I don't think we need it. */

/* 20. What is the state of setting useDefault in the various group */
/*     read/parse functions.  How many are actually ready for prime time? */
	
/* 22. radiation_parms reads randomParms from restart as a LITERAL.  Does */
/*     anyone have different/better ideas? */

/* 23. After talking to Jim some decisions are starting to clear up: */
/*     - For the single stream randoms Jim wants to get the same stream */
/*       of numbers on each task---however icky that may be.  Of course  */
/*       this begs the question which task's state do we save?  Unless */
/*       each task has called the generator the same number of times the */
/*       states will be different. */
/*     - Jim favors, and I concur, the notion of having only one scalar */
/*       state per generator instance.  All clients would share that */
/*       state.  If someone really wants to have their very own generator */
/*       they can ask for it in object.data */
/*     - We want to migrate to a design where RANDOM is directly */
/*       responsible for writing restart data rather than just serving */
/*       requests of others. */



/* Local Variables: */
/* tab-width: 3 */
/* End: */
