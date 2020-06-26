#include <assert.h>
#include <string.h>
#include "ddc.h"
#include "pio.h"
#include "units.h"
#include "io.h"
#include "mpiUtils.h"
#include "mpiTypes.h"
#include "ddcMalloc.h"
#include "domain.h"
#include "heap.h"
#include "recordDecorator.h"
extern FILE* ddcfile;

static char* appendTokens(char* s1, char* s2)
{
	if (s1 == NULL)
		s1 = strdup("");
	
	s1 = ddcRealloc(s1, strlen(s1) + strlen(s2) + 2);
	if (strlen(s1) > 0)
		strcat(s1, " ");
	strcat(s1, s2);
	return s1;
}


void ddc_writePXYZ(DDC *ddc, SIMULATE* simulate, PFILE*file)
{
   pio_long64 size = ddc->domains.size; 
   unsigned nremote = ddc->number_remote;
   unsigned nlocal = ddc->number_local;

   unsigned lrec = 96;
   unsigned nfields = 6;
   char* fieldTypes = NULL;
   char* fieldNames = NULL;
   char* fieldUnits = NULL;
	fieldTypes = appendTokens(fieldTypes, "u f f f u u");
	fieldNames = appendTokens(fieldNames, "id rx ry rz nlocal nremote");
	fieldUnits = appendTokens(fieldUnits, "1 Angstrom Angstrom Angstrom 1 1");
   char fmt[] = "%6u %20.13f %20.13f %20.13f %8d %8d";

	for (unsigned ii=0; ii<ddc->nPxyzDecorator; ++ii)
	{
		RECORD_DECORATOR* d = ddc->pxyzDecorator[ii];
		lrec      += d->lrec;
		nfields   += d->nfields;
		fieldTypes = appendTokens(fieldTypes, d->fieldTypes);
		fieldNames = appendTokens(fieldNames, d->fieldNames);
		fieldUnits = appendTokens(fieldUnits, d->fieldUnits);
	}
	   
   PioSet(file, "datatype", FIXRECORDASCII);
   PioSet(file, "recordLength", lrec);
   PioSet(file, "numberRecords", size);
   PioSet(file, "checksum", PIO_NONE);
   PioSet(file, "nfields", nfields);
   PioSet(file, "field_names", fieldNames);
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "field_units", fieldUnits);
   if (getRank(0) == 0)	write_fileheader(file, simulate, "pxyz");
	

   unsigned id = ddc->domains.local_id; 
   double lengthConvert = units_convert(1, "length_internal", "Angstrom");
   double x=ddc->domains.domains[id].center.x * lengthConvert;
   double y=ddc->domains.domains[id].center.y * lengthConvert;
   double z=ddc->domains.domains[id].center.z * lengthConvert;
   
   char buf[lrec];
   snprintf(buf, lrec, fmt, id, x, y, z, nlocal, nremote);
	for (unsigned ii=0; ii<ddc->nPxyzDecorator; ++ii)
	{
		strcat(buf, " ");
		ddc->pxyzDecorator[ii]->print(buf+strlen(buf), 0);
	}

	for (unsigned ii=strlen(buf); ii<lrec-1; ++ii)
      buf[ii] = ' ';
   buf[lrec-1] = '\n';
   Pwrite(buf, lrec, 1, file);
	ddcFree(fieldUnits);
	ddcFree(fieldNames);
	ddcFree(fieldTypes);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
