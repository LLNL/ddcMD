#ifndef THERMALIZE_H
#define THERMALIZE_H

#include "system.h"

enum THERMALIZE_ENUM 
{
	THERMALIZE_GLOBAL, THERMALIZE_BY_SPECIES, THERMALIZE_BY_GROUP,
	THERMALIZE_BOLTZMANN, THERMALIZE_JUTTNER, THERMALIZE_RESCALE
};

// To Do:
//
// 1.  The keepVcm option only works with the RESCALE method.

/**
 *
 *   # Temperatures are expected to be in internal units.
 */
typedef struct thermalize_parms_st
{
   enum THERMALIZE_ENUM selectionMode;
	enum THERMALIZE_ENUM method;
   double*  temperature;
   LONG64   seed;
	unsigned randomize;
   unsigned keepVcm;
} THERMALIZE_PARMS;

/** Construct a THERMALIZE_PARMS object and initialize all parms values
 *  with defaults.  The caller is responsible to call thermalize_destroy
 *  and to free the returned pointer. */
THERMALIZE_PARMS* thermalize_init(void);
void thermalize_destroy(THERMALIZE_PARMS* thisParms);


void thermalize(SYSTEM* system, THERMALIZE_PARMS* parms);

/** This is a convenience function that you can use to set the
 * temperature parameters in a THERMALIZE_PARMS object without worrying
 * about the internal details.  Just pass in an array of temperatures
 * and corresponding group or species names.  If the mode is
 * THERMALIZE_GLOBAL then the names are ignored (may be NULL) and only a
 * single temperature is used.  */
int thermalize_setTemperature(
	THERMALIZE_PARMS* thisParms, 
	enum THERMALIZE_ENUM mode, double* temperature,
	char** names, unsigned nNames);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
