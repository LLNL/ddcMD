#include "eq.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "ddcMalloc.h"
#include "units.h"
double eqconstant(double t, EQTARGET*eq)
{
	eq->value = eq->v0;
	return eq->v0;
}
double eqconstantIntegral(double t1, double t2, EQTARGET*eq)
{
	return eq->v0*(t2-t1);
}

double eqstep(double t, EQTARGET*eq)
{
	if (t < eq->t0) return eq->value = eq->v0;
	eq->value = eq->v1;
	return eq->v1;
}
double eqstepIntegral(double t1, double t2, EQTARGET*eq)
{
	double Integral1,Integral2; 
	if (t1 < eq->t0) Integral1 = eq->v0*t1; else Integral1 = eq->v1*t1; 
	if (t2 < eq->t0) Integral2 = eq->v0*t2; else Integral2 = eq->v1*t2; 
	return Integral2-Integral1;
}

double eqramp(double t, EQTARGET*eq)
{

	if (t < eq->t0) return eq->value = eq->v0;
	if (t > eq->t0 + eq->tau) return eq->value = eq->v1;
	return eq->value = eq->v0 + (eq->v1 - eq->v0)*(t - eq->t0)/eq->tau;
}
double eqrampIntegral(double t1, double t2, EQTARGET*eq)
{
	double Integral1,Integral2; 
	if (t1 < eq->t0) Integral1 = eq->v0*t1; 
	else if (t1> eq->t0+eq->tau) Integral1 = eq->v1*t1 ; 
	else  Integral1 = eq->v0*t1 + 0.5*(eq->v1 - eq->v0)*(t1 - eq->t0)*(t1-eq->t0)/eq->tau;

	if (t2 < eq->t0) Integral2 = eq->v0*t2; 
	else if (t2> eq->t0+eq->tau) Integral2 = eq->v1*t2 ; 
	else  Integral2 = eq->v0*t2 + 0.5*(eq->v1 - eq->v0)*(t2 - eq->t0)*(t2-eq->t0)/eq->tau;

	return Integral2-Integral1;

}

double eqexponential(double t, EQTARGET*eq)
{
	double f;
	f = exp((eq->t0 - t)/eq->tau);
	if (t < eq->t0) return eq->value = eq->v0;
	return eq->value = eq->v0*f + eq->v1*(1.0 - f);
}
double eqexponentialIntegral(double t1, double t2, EQTARGET*eq)
{
	double f;
	double Integral1,Integral2; 
	if (t1 < eq->t0) Integral1 = eq->v0*t1; 
	else 
	{
		f = exp((eq->t0 - t1)/eq->tau);
		Integral1 = -eq->tau*(eq->v0*f + eq->v1*(1.0 - f));
	}
	if (t2 < eq->t0) Integral2 = eq->v0*t2; 
	else 
	{
		f = exp((eq->t0 - t2)/eq->tau);
		Integral2 = -eq->tau*(eq->v0*f + eq->v1*(1.0 - f));
	}
	return Integral2-Integral1;

}
double eqcos(double t, EQTARGET*eq)
{
	double f;
	if (t < eq->t0) return eq->value = eq->v0;
	f = eq->value=0.5*( (eq->v0+eq->v1) + (eq->v0-eq->v1)*cos(2.0*M_PI*(t-eq->t0)/eq->tau) );
	return f ;
}
double eqcosIntegral(double t1, double t2, EQTARGET*eq)
{
	double Integral1,Integral2; 
	if (t1 < eq->t0) Integral1 = eq->v0*t1; 
	else 
	{
		Integral1 = eq->value=0.5*( (eq->v0+eq->v1)*t1 + (eq->tau/(2.0*M_PI)) * (eq->v0-eq->v1)*sin(2.0*M_PI*(t1-eq->t0)/eq->tau) );
	}
	if (t2 < eq->t0) Integral2 = eq->v0*t2; 
	else 
	{
		Integral2 = eq->value=0.5*( (eq->v0+eq->v1)*t2 + (eq->tau/(2.0*M_PI)) * (eq->v0-eq->v1)*sin(2.0*M_PI*(t2-eq->t0)/eq->tau) );
	}
	return Integral2-Integral1 ;
}

EQTARGET *eq_parse(char *string,char *return_units,char *arg_units)
{
	char *tok, *eptr;
	double v0,v1,t0,tau; 
	EQTARGET *eqtarget;
	eqtarget = ddcMalloc(sizeof(EQTARGET));
	eqtarget->function = NULL;

	prune_spaces(string); 
	zapChar('"', string);
	tok = strtok(string, " (");
	if (strcasecmp(tok, "RAMP") == 0){eqtarget->function = eqramp; eqtarget->integral = eqrampIntegral;  } 
	if (strcasecmp(tok, "STEP") == 0){eqtarget->function = eqstep; eqtarget->integral = eqstepIntegral; }
	if (strcasecmp(tok, "EXP") == 0) {eqtarget->function = eqexponential; eqtarget->integral = eqexponentialIntegral; }
	if (strcasecmp(tok, "COS") == 0) {eqtarget->function = eqcos; eqtarget->integral = eqcosIntegral; } 
	if (eqtarget->function != NULL)
	{
		tok = strtok(NULL, " ,");
		v0 = strtod(tok, &eptr); 
		trim(eptr); 
		if (*eptr != (char)0) v0 = units_convert(v0,eptr,NULL); 
		else v0 = units_convert(v0,return_units,NULL); 
		tok = strtok(NULL, " ,");
		v1 = strtod(tok, &eptr); 
		trim(eptr); 
		if (*eptr != (char)0) v1 = units_convert(v1,eptr,NULL); else v1 = units_convert(v1,return_units,NULL); 
		tok = strtok(NULL, " ,");
		t0 = strtod(tok, &eptr); trim(eptr); 
      if (*eptr != (char)0) t0 = units_convert(t0,eptr,NULL); else t0 = units_convert(t0,arg_units,NULL); 
		tok = strtok(NULL, " )");
		tau = strtod(tok, &eptr); trim(eptr); 
      if (*eptr != (char)0) tau = units_convert(tau,eptr,NULL); else tau = units_convert(tau,arg_units,NULL); 

		tok = strtok(NULL, " ,");
		eqtarget->v0 = v0;
		eqtarget->v1 = v1;
		eqtarget->t0 = t0;
		eqtarget->tau = tau;
		eqtarget->value = 0.0;
		return eqtarget;
	}
	eqtarget->function = eqconstant;
	eqtarget->integral = eqconstantIntegral;
	v0 = strtod(tok, &eptr); trim(eptr); if (*eptr != (char)0) v0 = units_convert(v0,eptr,NULL); else v0 = units_convert(v0,return_units,NULL); 
	eqtarget->value = eqtarget->v0 = v0;
	return eqtarget;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
