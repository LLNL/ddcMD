#ifndef PRINTINFO_H
#define PRINTINFO_H
#include "energyInfo.h"
#include "simulate.h"
enum PRINTINFO_ENUM { PLENGTH,PMASS,PTIME,PCURRENT,PTEMPERATURE,PAMOUNT,PLUMINOUS_INTENSITY,PENERGY,PPRESSURE,PVOLUME,PFORCE,PENERGYFLUX,NUNITS}; 
typedef struct units_st { char *name; char *value; double convert;} UNITS; 
typedef struct printinfo_st
{
	char *name;		/* name of the system */
	char *objclass;
	char *value;
	char *type;		/* type label */
	void  *parent; 
	int header;
	int printStress; // write stress.data?
	int printHmat;   // write hmatrix.data?
	int printGroup;   // write group info?
	int printGraphs;   // write graphs?
	int printMolecularPressure;   // write graphs?
	UNITS units[NUNITS]; 
} PRINTINFO;

PRINTINFO *printinfo_init(void *parent,char *name);
void printinfo_close(SIMULATE *simulate) ;
void printinfo(SIMULATE *simulate, ETYPE *energyInfo);
void printinfoAll(SIMULATE *simulate, ETYPE *energyInfo);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
