#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "error.h"
#include "three_algebra.h"
#include "object.h"
#include "pio.h"
#include "io.h"
#include "ioUtils.h"
#include "units.h"
#include "crc32.h"
#include "ddcMalloc.h"
#include "mpiAlgorithm.h"

int getRank(int);


typedef struct stressWrite_parms_st
{
   int    _nFiles;
   char*  _stressUnits;
} STRESS_WRITE_PARMS;


STRESS_WRITE_PARMS* stressWrite_parms(ANALYSIS *analysis)
{
   STRESS_WRITE_PARMS *parms = ddcMalloc(sizeof(STRESS_WRITE_PARMS)); 
	
   OBJECT* obj = (OBJECT*) analysis;

   object_get(obj, "nfiles", &parms->_nFiles, INT, 1, "0");
   object_get(obj, "units",  &parms->_stressUnits, LITERAL, 1, "eV/Angstrom^3");
   
   zapChar('"', parms->_stressUnits);
   if (units_check("pressure_internal", parms->_stressUnits) != 0)
      error_action(
	 "unit specified for stressWrite analysis must have units of pressure.",
	 ERROR_IN("stressWrite_parms", ABORT));
	      
   return parms;
}

void stressWrite_output(ANALYSIS* analysis)
{

   STRESS_WRITE_PARMS* parms = (STRESS_WRITE_PARMS *)analysis->parms; 
   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state; 

   gid_type* label = state->label;
   unsigned nlocal = state->nlocal;
   gid_type nglobal = sys->nglobal;

   THREE_SMATRIX* stress = state->sion;
   
   char fieldTypes[1024];
   unsigned gidFieldSize = bFieldSize(mpiMaxVal(label, nlocal,COMM_LOCAL));
   snprintf(fieldTypes, 1023, "u4 b%1u f4 f4 f4 f4 f4 f4", gidFieldSize);
   char fieldUnits[1024];
   snprintf(fieldUnits, 1023, "1 1 %s %s %s %s %s %s",
	    parms->_stressUnits, parms->_stressUnits, parms->_stressUnits,
	    parms->_stressUnits, parms->_stressUnits, parms->_stressUnits);
   unsigned offset[8];
   unsigned nParse = parseFieldSizes(offset, fieldTypes);
   assert(nParse == 8);
   
   CreateSnapshotdir(simulate, NULL);
   char filename[1024];
   snprintf(filename, 1023,"%s/stress", simulate->snapshotdir);
   PFILE* file = Popen(filename, "w", COMM_LOCAL);
   int lrec = gidFieldSize + 4 + 6*4; 
   unsigned char* line = (unsigned char*) ddcMalloc(lrec*sizeof(unsigned char));
   PioSet(file, "recordLength", lrec);
   PioSet(file, "datatype", FIXRECORDBINARY);
   PioSet(file, "numberRecords", nglobal);
   PioSet(file, "checksum", CRC32);
   PioSet(file, "nfields", 8);
   PioSet(file, "field_names", "checksum id sxx syy szz sxy sxz syz" );
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "field_units", fieldUnits);
   if (parms->_nFiles > 0)
      PioSet(file, "ngroup", parms->_nFiles);
   if (getRank(0) == 0) write_fileheader(file, simulate, "stress" );

   double pressureConvert = units_convert(1.0, "pressure_internal", parms->_stressUnits); 

   for (unsigned ii=0; ii<sys->nlocal; ++ii)
   {
      float f4;
      unsigned i4;
      unsigned n=1;
      bFieldPack(line+offset[n++], gidFieldSize, label[ii]);
      f4 = stress[ii].xx * pressureConvert; copyBytes(line+offset[n++], &f4, 4);
      f4 = stress[ii].yy * pressureConvert; copyBytes(line+offset[n++], &f4, 4);
      f4 = stress[ii].zz * pressureConvert; copyBytes(line+offset[n++], &f4, 4);
      f4 = stress[ii].xy * pressureConvert; copyBytes(line+offset[n++], &f4, 4);
      f4 = stress[ii].xz * pressureConvert; copyBytes(line+offset[n++], &f4, 4);
      f4 = stress[ii].yz * pressureConvert; copyBytes(line+offset[n++], &f4, 4);

      i4 =  checksum_crc32_table(line+offset[1],lrec-offset[1]);
      copyBytes(line, &i4, 4);
      Pwrite(line, lrec, 1, (PFILE *)file);	   
       
   }
   ddcFree(line);
   Pclose(file);
}

void stressWrite_eval(ANALYSIS* analysis)
{
   // empty
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
