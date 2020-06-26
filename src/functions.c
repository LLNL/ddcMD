#include  "functions.h"
#include  <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include <math.h>
#include "object.h"
#include "error.h"
#include "three_algebra.h"
#include "ddcMalloc.h"
double linear(LINEAR_PARMS*parms, double x, double *dfdx)
{
	*dfdx = parms->m;
	return parms->m*x + parms->b;
}

double stepdown(STEPDOWN_PARMS*parms, double x, double *dfdx)
{
	double x1, x2, f, dxi;
	*dfdx = 0.0;
	x1 = parms->x1;
	x2 = parms->x2;
	if (x <= x1) return 1.0;
	if (x >= x2) return 0.0;
	dxi = 1.0/(x2 - x1);
	f = 0.5 + 0.5*cos(M_PI*(x - x1)*dxi);
	*dfdx = -0.5*dxi*M_PI*sin(M_PI*(x - x1)*dxi);
	return f;
}

FUNCTION *function_init(void *parent, char *name)
{
	FUNCTION *function;
	function = (FUNCTION *) object_initialize(name, "FUNCTION", sizeof(FUNCTION));
	function->parent = parent; 
	object_get((OBJECT *) function, "type", &function->type, STRING, 1, "LINEAR");
	if (strcmp(function->type, "LINEAR") == 0)
	{
		LINEAR_PARMS *parms;
		function->itype = LINEAR;
		function->parms = parms = ddcMalloc(sizeof(LINEAR_PARMS));
		function->f = (double (*)())linear;
		object_get((OBJECT *) function, "m", &parms->m, DOUBLE, 1, "0.0");
		object_get((OBJECT *) function, "b", &parms->b, DOUBLE, 1, "0.0");
	}
	if (strcmp(function->type, "STEPDOWN") == 0)
	{
		STEPDOWN_PARMS *parms;
		function->itype = STEPDOWN;
		function->parms = parms = ddcMalloc(sizeof(STEPDOWN_PARMS));
		function->f = (double (*)())stepdown;
		object_get((OBJECT *) function, "x1", &parms->x1, DOUBLE, 1, "0.0");
		object_get((OBJECT *) function, "x2", &parms->x2, DOUBLE, 1, "0.0");
	}
	return function;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
