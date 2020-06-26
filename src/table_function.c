#include <math.h>
#include <stdlib.h>
#include "system.h"
#include "neighbor.h"
#include "object.h"
#include "units.h"
enum TABLE_ENUM { ARBITRARY_INTERVALS,UNIFORM_INTERVALS};
typedef struct table_element_struct 
{ 
	double x; double *a;
} TABLE_ELEMENT;
typedef struct table_struct 
{
	unsigned nspecies; 
	char **species_list; 
	RCUT_TYPE *rcut;
	double rmax;
	double(*fcn)(void *,double, double *); 
	int interval; 
	int number_intervals; 
	int number_terms; 
	double Rmin,Rmax;
	double interval_width; 
	TABLE_ELEMENT *element; 
	double *coeff;
	char *filename; 
}  TABLE_PARMS; 
double table_function_uniform(TABLE_PARMS *table,double r,double *dvdr) ;
TABLE_PARMS *table_parms(POTENTIAL *potential)
{
	char line[1024],*energyUnits,*lengthUnits; 
	TABLE_PARMS *table ;
	table =malloc(sizeof(TABLE_PARMS));
	SYSTEM *sys = (SYSTEM *)potential->parent;
	if (sys!= NULL) table->nspecies=sys->nspecies;  else table->nspecies = 1; 
	table->rcut  = malloc(sizeof(RCUT_TYPE)*1+(table->nspecies*table->nspecies)); 
	table->fcn=(double(*)(void*,double,double*))table_function_uniform;
	object_get((OBJECT *)potential,"number_intervals",&table->number_intervals,INT,1,"1"); 
	object_get((OBJECT *)potential,"number_terms",&table->number_terms,INT,1,"1"); 
	object_get((OBJECT *)potential,"filename",&table->filename,STRING,1,"table.data"); 
	object_get((OBJECT *)potential,"table_energyUnits",&energyUnits,STRING,1,"energy"); 
	object_get((OBJECT *)potential,"table_lengthUnits",&lengthUnits,STRING,1,"l"); 
	object_get((OBJECT *)potential, "Rmax", &table->Rmax, WITH_UNITS, 1, "0.0", "l", NULL);
	table->rmax=table->rcut[0].value = table->Rmax; 
	FILE *file  = fopen(table->filename,"r"); 
	table->element = malloc(sizeof(TABLE_ELEMENT)*table->number_intervals);
	table->coeff = malloc(sizeof(double)*table->number_intervals*table->number_terms);
	double energy_convert = units_convert(1.0,energyUnits,NULL); 
	double length_convert = units_convert(1.0,lengthUnits,NULL); 
	for (int i=0;i<table->number_intervals;i++)
	{
		char *ptr; 
		double convert = energy_convert; 
		table->element[i].a= table->coeff+i*table->number_terms; 
		fgets(line,1023,file); 
		table->element[i].x = strtod(line,&ptr);
		for (int j=0;j<table->number_terms;j++) 
		{
			table->element[i].a[j] = strtod(ptr,&ptr)*convert; 
			convert /= length_convert; 
		}
	}
	table->Rmin = table->element[0].x; 
	double s1=0.0,s2=0.0; 
	for (int i=0;i<table->number_intervals-1;i++)
	{
		double h = table->element[i+1].x-table->element[i].x; 
		s1 += h; 
		s2 += h*h; 
	}
	s1/= (table->number_intervals-1); 
	s2/= (table->number_intervals-1); 
	if ( fabs(1.0-s1*s1/s2) < 1e-12  ) 
	{
		table->interval = UNIFORM_INTERVALS; 
		table->interval_width= s1; 
	}
	else 
	{
		table->interval= ARBITRARY_INTERVALS; 
	}
	return table; 
}
double table_function_uniform(TABLE_PARMS *table,double r,double *dvdr) 
{
	if (r > table->Rmax) { *dvdr =0.0; return 0.0; }
	int i =  (r-table->Rmin)/table->interval_width; 
	TABLE_ELEMENT *element = table->element+i; 
	double d=0.0,v =0.0; 
	double x = r - element->x; 
	for (i=table->number_terms-1;i>0;i--)
	{
		v = element->a[i] + x*v; 
		d = i*element->a[i] + x*d; 
	}
	v = element->a[0] + x*v; 
	*dvdr = d; 
	return v; 
}
/*
#include "codata.h"
int main()
{
	int i; 
	double v1;
	units_internal(a0_MKS,Rinfhc_MKS*1e-30/(a0_MKS*a0_MKS),1e-15,e_MKS/1e-15,Rinfhc_eV/kB_eV,1.0,1.0);
	units_external(1e-10,u_MKS,1e-15,e_MKS/1e-15,1.0,1.0,1.0);
	object_compile();
	TABLE_PARMS *table; 
	POTENTIAL *obj = object_initialize("pair", "POTENTIAL", sizeof(TABLE_PARMS));
	obj->parent = NULL; 
	table = table_parms(obj); 
	//for (double r=0.0;r<16.0;r+=0.001)
	double r,v; 
	double convert = units_convert(1.0,"l",NULL); 
	double convert_energy = units_convert(1.0,NULL,"Hartree"); 
	double convert_pressure = units_convert(1.0,"eV/Ang^3","GPa"); 
	scanf("%lf %le",&r,&v1);
	double	rho = 1.0/(128.450190688425)*convert*convert*convert; 
	while (!feof(stdin))
	{
		double dvdr,v; 
		v=convert_energy*table_function_uniform(table,r,&dvdr);
		dvdr *= (convert_energy*convert);
		printf("%20.14f %24.14e %24.14e %24.14e %24.14e\n",r,v,v1,fabs(v-v1),(v-v1)/(fabs(v1)+1e-8)); 
		scanf("%lf %le",&r,&v1);
	}
//	printf("1 ev/Ang^3=%e %e\n",convert_pressure,rho*1000/10604.5*convert_pressure); 
}
*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
