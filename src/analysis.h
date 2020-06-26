#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "gid.h"

enum analysis_enum {
   NO_ANALYSIS, COARSEGRAIN, SUBSET_WRITE, VELOCITYAUTOCORRELATION,
   PAIRCORRELATION, STRESS_WRITE, STATIC_STRUCTURE_FACTOR,
   DYNAMIC_STRUCTURE_FACTOR, VCM_WRITE, PLASMAWAKE, FORCEAVERAGE,
   CENTRO_SYMMETRY,BOUND_STATE, ACKLAND_JONES, VON_MISES, QUATERNION,
   NEAR_MISS, ZDENSITY,KINETICENERGYDISTN, EWALD_TRANSPORT_COEFF, DATASUBSET,
   PROJECTILE_RV,PARTICLE_HISTOGRAM,LOCALYUKAWA_ANALYSIS, EWALD_ANALYSIS, 
   EVENTS, CDB_ANALYSIS, CHOL_ANALYSIS, PAIRANALYSIS
};

/** Implementation notes for analysis functions:
 *
 *  Every analysis object must implement three functions: a parms
 *  function to get data from the object database, an eval function to
 *  perform the analysis, and a write function to write the analyzed
 *  data to disk.  Analysis objects may also optionally implement a
 *  clear function that discards all previously accumulated data.  In
 *  cases where there is no method to accumulate and average data over
 *  time the eval and write functions can be combined, however both
 *  functions must still be provided (one can be left empty).
 *
 *  The eval function will be called every eval_rate time steps.  It
 *  will also be called once when the -RW option is specified on the
 *  command line.
 *
 *  The write function is called every outputrate time steps.  It
 *  will also be called once when the -RW option is specified on the
 *  command line.
 *
 *  The clear function is called before entering the main time step
 *  loop.
 *
 *  Implementers must choose whether to check the eval and output rates
 *  in their eval and write implementations.  It may be desirable to
 *  not check so that the -RW switch can always result in analysis
 *  output.
 *
 *  Note that the eval and output rate variables are in the analysis
 *  base class are are scanned from the object database by the
 *  analysis_init code.  Subclass authors need not worry about scanning
 *  these parameters.
 */


typedef struct analysis_st
{
	char *name;		/* analysis name */
	char *objclass;
	char *value;
	char *type;		/* model */
	void *parent; 
	int  itype;	/* integer label for type */
	void (*eval) (struct analysis_st  *a);	/* pointer to analysis function */
	void (*output) (struct analysis_st  *a);	/* pointer to analysis output function */
	void (*clear) (struct analysis_st  *a);	/* pointer to analysis clear */
	void (*startup) (struct analysis_st  *a);	/* pointer to analysis clear */
	int eval_rate,outputrate ; 
	void *parms;		/* pointer to  parameters for analysis function */
} ANALYSIS;



#endif // #ifndef ANALYSIS_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
