#ifndef EQ_H
#define EQ_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct eqtarget_st
{
	double (*function) (double t, struct eqtarget_st *);
	double (*integral) (double t1, double t2, struct eqtarget_st *);
	double v0, v1, t0, tau;
	double value;
} EQTARGET;
EQTARGET *eq_parse(char *string, char *, char *);
double eqconstant(double t, EQTARGET*eq);
double eqconstantIntegral(double t1, double t2, EQTARGET*eq);
double eqstep(double t, EQTARGET*eq);
double eqstepIntegral(double t1, double t2, EQTARGET*eq);
double eqramp(double t, EQTARGET*eq);
double eqrampIntegral(double t1, double t2, EQTARGET*eq);
double eqexponential(double t, EQTARGET*eq);
double eqexponentialIntegral(double t1, double t2, EQTARGET*eq);
double eqcos(double t, EQTARGET*eq);
double eqcosIntegral(double t1, double t2, EQTARGET*eq);

#ifdef __cplusplus
}
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
