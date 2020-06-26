#include "nearPair.h"
#include "error.h"
#include "expandbuffer.h"

static NEARPAIR_PARMS nearPairParms={0,0,0,0,0};  

NEARPAIR_PARMS getNearPairParms()
{
    return nearPairParms; 
}
NEARPAIR *getNearPairs()
{
    return nearPairParms.pairs; 
}
int getNearPairNumber()
{
    return nearPairParms.number; 
}
void  nearPairSetAlpha(double alpha) 
{
	nearPairParms.alpha=alpha; 
}

void setNearPairRmax(double rmax)
{
	if (rmax > nearPairParms.rmax) nearPairParms.rmax=rmax; 
}
double getNearPairRmax()
{
	return nearPairParms.rmax; 
}
void addPair2List(int i, int j, double r) 
{
	
   nearPairParms.pairs  = (NEARPAIR *) ExpandBuffers((void *)nearPairParms.pairs, sizeof(NEARPAIR), nearPairParms.number+1, 64, LOCATION("addPair2List"),"nearPairParms->pair");
   nearPairParms.pairs[nearPairParms.number].i=i; 
	nearPairParms.pairs[nearPairParms.number].j=j;
	nearPairParms.pairs[nearPairParms.number].r=r; 
	nearPairParms.number++; 
}
