#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "gid.h"

enum transform_enum {NO_TRANSFORM, VELOCITY_TRANSFORM, SELECT_TRANSFORM,
		     REPLICATE_TRANSFORM, THERMALIZE_TRANSFORM, TRANSECTMORPH_TRANSFORM,
           ALCHEMY_TRANSFORM, BOX_TRANSFORM,
           LINEARISOTROPICV_TRANSFORM, GID_SHUFFLE_TRANSFORM,
           PROJECTILE_TRANSFORM, SHOCK_TRANSFORM, APPEND_TRANSFORM, ASSIGNGROUPS_TRANSFORM, IMPACT_TRANSFORM, CUSTOM_TRANSFORM}; 

enum transformRate_enum { atMaxLoop=-5, atStartThenExit = -4, atStart = -3, atCheckpoint = -2, atFinish = -1, never =0};

typedef struct transform_st
{
   char *name;		/* transform name */
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   int  itype;	       /* integer label for type */
   void (*function) (struct transform_st  *a);	/* ptr to transform function */
   void *parms;		/* ptr to  parameters for transform function */
   char* writedir;
   SIGNED64  loop;
   double time;
   int rate; 
   int writeNeeded;
} TRANSFORM;


/** Default values:
 *  The default time and loop are the time and loop count read from the
 *  restart file.  The default checkpoint name is transform.n where n is
 *  the 8 loop count written as an 8 digit number.
 */
TRANSFORM* transform_init(void* parent, char* name);

/** The applyTranforms functions iterates the array of transforms and
 * applies them sequentially in the order they are specified in the
 * SIMULATE object.  After all transforms are applied a checkpoint file
 * is written with the name, time, and loop that are specified in the
 * last transform.  EXCEPTION:  Some transforms (such as replicate)
 * perform their own save operations.  Such transforms will set
 * transform->writeNeeded=0 to inhibit the writing in the
 * applyTransforms function.  Such transforms *must* be the last in a
 * sequence or any subsequent transforms will be ignored.
 **/
void atCheckpointTransforms(int nTransform, TRANSFORM** transform);
void atStartThenExitTransforms(int nTransform, TRANSFORM** transform);
void atRateTransforms(int nTransforms, TRANSFORM** transform);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
