#ifndef EAM_H
#define EAM_H
enum SERIES_ORDER  {  NO_INTERLACE, PARTIAL_INTERLACE, FULL_INTERLACE } ;
typedef struct ep_st { double e,p; } EP;


/** If this were C++ this would be a base class and the PARMS classes
 *  for other implementations of the EAM potential such as EAM_OPT or
 *  EAM_ONEPASS would be polymorphic subclasses.  Since this is C, those
 *  other PARMS classes just copy the members of this class in identical
 *  order so that we can pass and EAM_FOO_PARMS* object to the various
 *  eam_FORM_parms functions regardless of what FOO is.
 *
 *  DO NOT CHANGE THIS STRUCT IN ANY WAY UNLESS YOU FIND ALL OF THE
 *  POLYMORPHIC SUBCLASSES AND MAKE CORRESPONDING CHANGES.
 *
 *  This structure represents all of the function calls that are
 *  necessary to implement a two-pass evaluation of a general EAM type
 *  potential.  The various methods of iterating over neighbors can call
 *  these functions thus separating the neighbor finding algorithm from
 *  the potential evaluation.
 */
typedef struct eam_parms_str
{
	double rmax;
	int nspecies; 
	double (*embedding_function)(void *parms, double rho, double *dv_dr);
	EP (*pass1_function)(void  *parms,double r);
	EP (*pass2_function)(void  *parms,double r);
	void **embedding_parms; 
	void **pass_parms; 
	int cmax;
	enum SERIES_ORDER SeriesOrder; 
} EAM_PARMS;
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
