#ifndef NEARPAIR_H
#define NEARPAIR_H
#include "three_algebra.h"
#include "gid.h"

typedef struct nearpair_st
{
	gid_type iLabel,jLabel; 
	int i,j, turningflag;
	double r,tmin,rmin,vmin,emin;
	THREE_VECTOR Fr;
	THREE_VECTOR rDelta0, vDelta0,rDelta,vDelta,rCM,vCM;
} NEARPAIR; 

typedef struct nearpair_parms_st
{
	double smallBallRadius;
	double rmax;
	double alpha;
	int number;
	NEARPAIR *pairs;
} NEARPAIR_PARMS; 

NEARPAIR_PARMS getNearPairParms();
NEARPAIR *getNearPairs();
double getNearPairRmax();
void setNearPairRmax(double rmax);
int getNearPairNumber();
void  nearPairSetAlpha(double alpha) ;
void addPair2List(int i, int j, double r); 
#endif 

/* Local Variables: */
/* tab-width: 3 */
/* End: */
