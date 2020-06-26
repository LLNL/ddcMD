#include "printinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>		/*mallinfo call  MEMORY DEBUG */
#include "object.h"
#include "md.h"
#include "three_algebra.h"
#include "mgpt.h"
#include "ddc.h"
#include "utilities.h"
#include "simulate.h"
#include "energyInfo.h"
#include "ddcMalloc.h"
#include "units.h"
#include "mpiUtils.h"
#include "orderSH.h"
#include "nextfile.h"
#include "nptglf.h"
#include "format.h"
#include "system.h"

static FILE *datafile = NULL;
static FILE *stressfile = NULL;
static FILE *hmatfile = NULL;
static FILE **groupdatafile = NULL;
static void print_stress(SIMULATE*simulate, char *name, ETYPE*energyInfo, FILE*datafile,int  header);
static void print_hmat(SIMULATE* simulate, FILE* datafile, int  header);
static char *unitnames[] =         { "LENGTH","MASS","TIME","CURRENT","TEMPERATURE","AMOUNT","LUMINOUS_INTENSITY","ENERGY","PRESSURE","VOLUME","FORCE","ENERGYFLUX"}; 
static char *unitDefaultValue[] = { "Ang",   "amu", "fs",  "e/fs",   "K",          "1",     "1",                  "eV",    "GPa",     "Bohr^3","eV/Ang","ueV/Ang^2/fs"}; 
THREE_SMATRIX molecularPressure(SYSTEM *sys,THREE_SMATRIX virial,double T);
PRINTINFO *printinfo_init(void *parent,char *name)
{
	static PRINTINFO *printinfo;
	char *type;
	char string[64];
	sprintf(string,"%s PRINTINFO {}",name);

	
	object_compilestring(string ) ;
	printinfo = (PRINTINFO *) object_initialize(name, "PRINTINFO", sizeof(PRINTINFO));
	printinfo->parent = parent; 

	printinfo->header = 1;
	object_get((OBJECT *) printinfo, "type", &(type), STRING, 1, "NORMAL");
	object_get((OBJECT*) printinfo, "printStress",  &printinfo->printStress, INT, 1, "0");
	object_get((OBJECT*) printinfo, "printHmatrix", &printinfo->printHmat,   INT, 1, "0");
	object_get((OBJECT*) printinfo, "printGroup", &printinfo->printGroup,   INT, 1, "0");
	object_get((OBJECT*) printinfo, "printGraphs", &printinfo->printGraphs,   INT, 1, "0");
	object_get((OBJECT*) printinfo, "printMolecularPressure", &printinfo->printMolecularPressure,   INT, 1, "0");
	printinfo->name = ddcCalloc(strlen(name) + 1, sizeof(char));
	printinfo->type = ddcCalloc(strlen(type) + 1, sizeof(char));
	strcpy(printinfo->name, name);
	strcpy(printinfo->type, type);
	for (int i=0;i<NUNITS;i++) 
	{
		printinfo->units[i].name=unitnames[i];
		object_get((OBJECT*) printinfo,unitnames[i],&printinfo->units[i].value,STRING,1,unitDefaultValue[i]);
//		printinfo->units[i].value=NULL;
	}

/*
	printinfo->units[PLENGTH].value=strdup("Ang");
	printinfo->units[PTIME].value=strdup("fs");
	printinfo->units[PTEMPERATURE].value=strdup("eV/kB");
	printinfo->units[PENERGY].value=strdup("eV");
	printinfo->units[PPRESSURE].value=strdup("eV/Ang^3");
	printinfo->units[PVOLUME].value=strdup("Ang^3");
	printinfo->units[PFORCE].value=strdup("eV/Ang");
	printinfo->units[PENERGYFLUX].value=strdup("ueV/Ang^2/fs");
*/

	for (int i=0;i<NUNITS;i++) 
	{
		if (printinfo->units[i].value != NULL) 
		{
			printinfo->units[i].convert=units_convert(1.0,NULL,printinfo->units[i].value);
		}
	}
	return printinfo;
}
void printinfo_close(SIMULATE *simulate) 
{
	int i; 
	if (getRank(0) != 0 ) return ; 
	fclose(datafile); 
	return; 
	for (i=0;i<simulate->system->ngroup;i++) 
	{
		if (groupdatafile[i]!=NULL) fclose(groupdatafile[i]); 
	}
	ddcFree(groupdatafile); 
}
void printinfo(SIMULATE*simulate, ETYPE*energyInfo)
{
	UNITS *units=simulate->printinfo->units; 
	THREE_SMATRIX sion;
	SYSTEM *sys;
	if (getRank(0) != 0) return;
	sys = simulate->system;
	sion = energyInfo->sion;
	//double peq = 0.0;
	//if (simulate->Peq != NULL) peq = simulate->Peq->value;
	SMATSCALE(sion, units[PPRESSURE].convert);
	printf("Natoms = %"PRIu64"\n", sys->nglobal);
	printf("energyInfo->eion %f  (%s) energyInfo->pion %f (%s)\n", units[PENERGY].convert*energyInfo->eion/sys->nglobal, units[PENERGY].value, units[PPRESSURE].convert*energyInfo->pion, units[PPRESSURE].value);
	printf("stress tensor (%s)\n",units[PRESSURE].value);
	printf("%f %f %f \n", sion.xx, sion.xy, sion.xz);
	printf("\t %f %f \n", sion.yy, sion.yz);
	printf("\t \t %f \n", sion.zz);
	fflush(stdout);
}
static char *scat512(char *a,char*b)
{	
	static char s[512]; 
	snprintf(s,511,"%s(%s)",a,b); 
	return s; 
}
void printinfoA(SIMULATE*simulate, char *name, ETYPE*energyInfo, FILE*datafile,int  header)
{
	if (getRank(0) != 0) return;
	UNITS *units=simulate->printinfo->units; 
	SYSTEM *sys = simulate->system;
	double Eref = 0.0;
	//if (sys->potential[0]->itype == MGPT) Eref = ((MGPT_PARMS *) sys->potential[0]->parms)->evol0;
	double eBath = energyInfo->eBath*(units[PENERGY].convert/sys->nglobal);
	SIGNED64 loop = simulate->loop;
	double time = units[PTIME].convert*simulate->time;
	double ekinetic = units[PENERGY].convert*(energyInfo->rk/sys->nglobal);
	double etot = units[PENERGY].convert*((energyInfo->eion + energyInfo->rk+energyInfo->eBath)/sys->nglobal);
	double epot = units[PENERGY].convert*(energyInfo->eion/sys->nglobal - Eref);
	double temperature = units[PTEMPERATURE].convert*energyInfo->temperature;
	double pressure = units[PPRESSURE].convert*energyInfo->pion;
	double voln = units[PVOLUME].convert*sys->box->volume/sys->nglobal;
	THREE_SMATRIX sion = energyInfo->sion; 
   THREE_MATRIX h0 = box_get_h(NULL); 
   MATSCALE(h0,units[PLENGTH].convert); 
	SMATSCALE(sion, units[PPRESSURE].convert);
   int i; 
	for (i = 0; i < sys->ngroup; i++) 
	{
		if (sys->group[i] && (sys->group[i]->itype == RADIATION || sys->group[i]->itype == IONIZATION || sys->group[i]->itype==TNBURN)) break; 
	}
	int writeEBath = ( i< sys->ngroup ); 
   int writeBoxDia   =1 ; 
	if (header)
	{
		char hstring[512];
		*hstring = (char)0; 
      char fmt[8]; sprintf(fmt,"%%-%ds",loopFormatSize()); 
		sprintf(hstring+strlen(hstring),fmt,"#loop"); 
		sprintf(hstring+strlen(hstring)," %16s",scat512("time", units[PTIME].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Etotal",    units[PENERGY].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Ekin",    units[PENERGY].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Epot",    units[PENERGY].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Temp",    units[PTEMPERATURE].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Press",    units[PPRESSURE].value));
		sprintf(hstring+strlen(hstring)," %18s", scat512("Volume",    units[PVOLUME].value));
      if (writeBoxDia) 
      {
		sprintf(hstring+strlen(hstring)," %15s", scat512("lx",    units[PLENGTH].value));
		sprintf(hstring+strlen(hstring)," %15s", scat512("ly",    units[PLENGTH].value));
		sprintf(hstring+strlen(hstring)," %15s", scat512("lz",    units[PLENGTH].value));
      }
		if (writeEBath) sprintf(hstring+strlen(hstring),"  %14s", scat512("EnergyBath",    units[PENERGY].value));
		for (i = 0; i < sys->npotential; i++)
		{
			if (sys->potential[i] && sys->potential[i]->itype == ORDERSH)
			{
				sprintf(hstring+strlen(hstring),"  %15s", "Phi");
				sprintf(hstring+strlen(hstring),"  %15s", scat512("EPhi",  units[PENERGY].value));
			}
		}
		if (simulate->integrator->itype == NPTGLF)
		{
				sprintf(hstring+strlen(hstring),"  %16s", scat512("P-Ptarget",  units[PPRESSURE].value));
				sprintf(hstring+strlen(hstring),"  %16s", "zeta");
		}
		printf("%s\n",hstring); 
		fprintf(datafile,"%s\n", hstring);
	}
	if (strcmp("System",name) ==0) 
	{
      printf(loopFormat(),loop);
		printf(" %16.6f %18.12f %18.12f %18.12f %18.8f %18.12f %18.12f ", time, etot, ekinetic, epot, temperature, pressure, voln);
  	if (writeBoxDia)	printf("%15.8f %15.8f %15.8f", h0.xx,h0.yy,h0.zz);
  	if (writeEBath)	printf("%15.12f ", eBath);
	}
   fprintf(datafile,loopFormat(),loop);
	fprintf(datafile, " %16.6f %18.12f %18.12f %18.12f %18.8f %18.12f %18.12f ", time, etot, ekinetic, epot, temperature, pressure, voln);
  	if (writeBoxDia)	fprintf(datafile,"%15.8f %15.8f %15.8f", h0.xx,h0.yy,h0.zz);
  	if (writeEBath)	fprintf(datafile,"%15.12f ", eBath);
	for (i = 0; i < sys->npotential; i++)
	{
		if (sys->potential[i] && sys->potential[i]->itype == ORDERSH)
		{
			ORDERSH_PARMS *parms;
			parms = (ORDERSH_PARMS *) (sys->potential[i]->parms);
			printf("%15.12f %15.12f ", parms->phi, parms->Energy*units[PENERGY].convert);
			fprintf(datafile, "%15.12f %15.12f ", parms->phi, parms->Energy*units[PENERGY].convert);
		}
	}
	if (simulate->integrator->itype == NPTGLF)
	{
		NPTGLF_PARMS *p;
		double zeta;
		double deltap;
		p = (NPTGLF_PARMS *) simulate->integrator->parms;
		//deltap = pressure - units[PPRESSURE].convert*simulate->Peq->value;
		deltap = pressure - units[PPRESSURE].convert*p->pressure;
		zeta = p->zeta;
		printf("%16.12f %16.12f ", deltap, zeta);
		fprintf(datafile, "%16.12f %16.12f ", deltap, zeta);
	}
	if (simulate->integrator->itype == NEXTFILE)
	{
		NEXTFILE_PARMS *p;
		p = (NEXTFILE_PARMS *) simulate->integrator->parms;
		printf("%s ", p->currentfile);
		fprintf(datafile, "%s ", p->currentfile);
	}
	if (strcmp("System",name) ==0) printf("\n");
	fprintf(datafile, "\n");
	fflush(stdout);
	fflush(datafile);
}
void convertToMolecularPressures(ETYPE *energyInfo)
{
    SYSTEM *sys=system_getSystem(NULL);
    THREE_SMATRIX pTensor = molecularPressure(sys,energyInfo->virial,energyInfo->temperature);
    energyInfo->pion = TRACE(pTensor)/3.0; 
    energyInfo->sion=pTensor;   
    SMATSCALE(energyInfo->sion,-1.0);
}
void printinfoAll(SIMULATE *simulate, ETYPE *energyInfo)
{
	if (getRank(0)  != 0 ) return;
	PRINTINFO* printinfo = simulate->printinfo;
	if (printinfo->header) 
	{
		datafile = fopen("data", "a");
		if (printinfo->printStress != 0) stressfile = fopen("stress.data", "a");
		if (printinfo->printHmat != 0)   hmatfile = fopen("hmatrix.data", "a");
	}
	SYSTEM* sys=simulate->system; 
   ETYPE energyInfoCpy=*energyInfo; 
   if (printinfo->printMolecularPressure) convertToMolecularPressures(&energyInfoCpy); 
	printinfoA(simulate, "System", &energyInfoCpy, datafile,printinfo->header);
	if (printinfo->printStress) print_stress(simulate, "System", energyInfo, stressfile,printinfo->header);
	if (printinfo->printHmat)   print_hmat(simulate, hmatfile, printinfo->header);
   if (printinfo->printGroup) 
   {
      if (sys->ngroup>1  )
      {
         if (printinfo->header) 
         {
            groupdatafile = ddcMalloc(sys->ngroup*sizeof(FILE *));
            for (int i=0;i<sys->ngroup;i++)
            {
               char filename[256];
               sprintf(filename,"%s.data",sys->group[i]->name);
               groupdatafile[i] = fopen(filename,"a"); 
            }
         }
         for (int i=0;i<sys->ngroup;i++)
         {
            if (sys->group[i]->energyInfo.number > 0.0) 
            {
               printinfoA(simulate, sys->group[i]->name,&sys->group[i]->energyInfo, groupdatafile[i],printinfo->header);
            }
         }
      }
   }
   printinfo->header=0; 
}
void print_stress(SIMULATE*simulate, char *name, ETYPE*energyInfo, FILE*file,int  header)
{
   UNITS *units=simulate->printinfo->units; 
   double time;
   SIGNED64 loop;
   if (getRank(0) != 0) return;
   loop = simulate->loop;
   time = units[PTIME].convert*simulate->time;
   if (header)
   {
      char hstring[512];
      *hstring = (char)0; 
      char fmt[8]; sprintf(fmt,"%%-%ds",loopFormatSize()); 
		sprintf(hstring+strlen(hstring),fmt,"#loop"); 
      sprintf(hstring+strlen(hstring)," %16s",scat512("time", units[PTIME].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_xx",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_yy",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_zz",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_xy",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_xz",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Sigma_yz",    units[PPRESSURE].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Jx",    units[PENERGYFLUX].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Jy",    units[PENERGYFLUX].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("Jz",    units[PENERGYFLUX].value));
      fprintf(file,"%s\n", hstring);
   }
   THREE_SMATRIX stress=energyInfo->sion; 
   THREE_VECTOR J=energyInfo->thermal_flux; 
   SMATSCALE(stress,units[PPRESSURE].convert); 
   VSCALE(J,units[PENERGYFLUX].convert); 

   fprintf(file, loopFormat(), loop); 
   fprintf(file, " %16.6f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f ", time, stress.xx,stress.yy,stress.zz,stress.xy,stress.xz,stress.yz);
   fprintf(file, "%16.12f %16.12f %16.12f\n", J.x,J.y,J.z);
   fflush(stdout);
   fflush(file);
}

void print_hmat(SIMULATE* simulate, FILE* file, int  header)
{
   UNITS *units=simulate->printinfo->units; 
   if (getRank(0) != 0) return;
   SIGNED64 loop = simulate->loop;
   double time = units[PTIME].convert*simulate->time;
   if (header)
   {
      char hstring[512];
      *hstring = (char)0; 
      char fmt[8]; sprintf(fmt,"%%-%ds",loopFormatSize()); 
      sprintf(hstring+strlen(hstring),fmt,"#loop"); 
      sprintf(hstring+strlen(hstring)," %16s", scat512("time", units[PTIME].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_xx", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_xy", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_xz", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_yx", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_yy", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_yz", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_zx", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_zy", units[PLENGTH].value));
      sprintf(hstring+strlen(hstring)," %16s", scat512("h_zz", units[PLENGTH].value));
      fprintf(file,"%s\n", hstring);
   }
   THREE_MATRIX h = simulate->system->box->h0;
   MATSCALE(h, units[PLENGTH].convert); 

   fprintf(file, loopFormat(), loop); 
   fprintf(file, " %16.6f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f\n",
         time, h.xx, h.xy, h.xz, h.yx, h.yy, h.yz, h.zx, h.zy, h.zz);
   fflush(file);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
