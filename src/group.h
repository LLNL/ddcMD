#ifndef GROUP_H
#define GROUP_H
#include <stdio.h>
#include "energyInfo.h"
#include "state.h"
#include "object.h"
#include "gid.h"

/**
 *  Group is an abstract base class for operations such as thermostats
 *  that modify particle dynamics outside of the usual F=ma.
 *
 *  Each derived group *must* define functions to resolve the following
 *  function pointers:
 *
 *    void     (*Update)        (struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*velocityUpdate)(int mode, int k, struct group_st *g, STATE *state, double time, double dt);
 *
 * Groups may also define the following optional functions:
 *    void (*defaultValue)      (struct group_st *g, gid_type label, int i);
 *    int  (*checkValue)         (struct group_st *g, int i);
 *    char*     (*parse)        (struct group_st *g, char *string, int i);
 *    char*    (*write)         (struct group_st *g, int i);
 *    unsigned (*bread)         (unsigned char* buf, struct group_st* g, int i);
 *    unsigned (*bwrite)        (unsigned char* buf, struct group_st *g, int i);
 *    void     (*write_dynamics)(struct group_st *g, FILE *file);
 *    void     (*start)         (struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*Update1)(struct group_st *g, int mode, void *state, double t, double dt);
 *    void     (*Update2)(struct group_st *g, int mode, void *state, double t, double dt);
 *
 *  By convention, each derived group also implements a parms function
 *  to scan the concrete class parameters from the input object and
 *  store them.
 *
 *  Update function:
 *  Called three times each timestep by the integrator.  The purpose is to
 *  update group parameters, private or static variables, etc. according
 *  to the state of the system.  Note that the timestep passed into the
 *  Update function is half of the actual timestep.  The mode argument
 *  will be set to either FRONT_TIMESTEP,BACK_TIMESTEP or BACK2_TIMESTEP
 *  depending on whether it is the first or second call (respectively) per
 *  timestep
 *
 *  velocityUpdate function:
 *  Called twice per timestep by the integrator to update the velocities
 *  of the particles.  Note that this function should add the quantity
 *  Newton's law quantity F*dt/m to each velocity (see
 *  free_velocityUpdate for an example) in addition to whatever
 *  modifications are required.
 *
 *  write and bwrite functions:
 *  Handles writing the per-atom data (such as the random number state
 *  for the langevin thermostat) into the ascii or binary (respectively)
 *  atoms file.  The write function returns the string to be written,
 *  the bwrite function populates the buffer provided in the argument
 *  list and returns the number of charaters written.
 *
 *  
 *  
 *  parse and bread functions
 *  Handles reading per-atom group data from the ascii and binary
 *  (respectively) atoms file.  In parse if read data is "bad" a flag
 *  useDefault is set to one.  (The bread function should do the same,
 *  but current implementations probably don't.)
 *
 *  defaultValue functions
 *  if useDefault flag is set (see above) the useDefault functions is called 
 *  by collection_init after read is completed. `
 *
 *  write_dynamics function
 *  Called by writeRestart to write group specific information into the
 *  restart file.
 *
 *  BUGS:
 *  There is currently no good method to check that the data in the
 *  atoms file corresponds to the actual group type.  This can cause all
 *  sorts of problems.
 */

enum GROUP_CLASS { ALL, FREE, RELATIVISTIC, FROZEN, NOSE, ANDERSON, LANGEVIN,
                   RELATIVISTICLANGEVIN, EXTFORCE, SHEAR, SHWALL, FIXEDVELOCITY,
                   RADIATION, IONIZATION, BERENDSEN, DOUBLE_MIRROR, UNIONGROUP,
                   PISTON, QUENCH, TNBURN, SHOCK};
enum GROUP_UPDATEMODE { START_GROUP, FRONT_TIMESTEP, BACK_TIMESTEP, BACK2_TIMESTEP} ;
enum GROUP_ENUM { GROUPTYPE, GROUPITYPE, GROUPLIST, NGROUPS, GROUPMAXNAMELENGTH, GROUPMAXWRITELENGTH, GROUPMAXBINARYWRITELENGTH};

typedef struct group_st
{
   char *name;		/* group name */
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   enum GROUP_CLASS itype; /* integer label for type */
   int index;
   void *parm;		/* pointer to  parameter */
   void *first, *last;
   gid_type  nMember;
   int useDefault; 
   ETYPE energyInfo; 
   double energyBath;
   char *field_names, *field_type, *field_units, *field_format; 
   void (*defaultValue) (struct group_st *g, gid_type label, int i);
   int (*checkValue) (struct group_st *g, int i);
   char *(*parse) (struct group_st *g, char *string, int i);
   char *(*write) (struct group_st *g, int i);
   unsigned (*bread) (unsigned char* buf, struct group_st* g, int i);
   unsigned (*bwrite) (unsigned char* buf, struct group_st *g, int i);
   void (*write_dynamics) (struct group_st *g, FILE *file);
   void (*start)(struct group_st *g, int mode, void *state, double t, double dt);
   void (*Update)(struct group_st *g, int mode, void *state, double t, double dt);
   void (*Update1)(struct group_st *g, int mode, void *state, double t, double dt);
   void (*Update2)(struct group_st *g, int mode, void *state, double t, double dt);
   void (*velocityUpdate)(int mode,int k, struct group_st *g, void  *state, double time, double dt);
   void (*velocityUpdateKernel)(int mode, struct group_st *g, void  *state, double time, double dt);
   double (*energyBathFunction)(struct group_st *group); 
} GROUP;
GROUP *group_init(void *parent,char *name);
void group_get(GROUP*group, int get, void **ptr);
int group_put(GROUP*group, int put, void **ptr);
GROUP *group_find(GROUP ** group, char *name);
GROUP *group_by_index(GROUP ** group, int index);
#endif 

/* Local Variables: */
/* tab-width: 3 */
/* End: */
