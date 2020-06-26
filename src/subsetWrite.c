#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "object.h"
#include "pio.h"
#include "io.h"
#include "crc32.h"
#include "ddcMalloc.h"
#include "units.h"
#include "mpiAlgorithm.h"
#include "pinfo.h"
#include "ioUtils.h"
#include "system.h"
#include "mpiUtils.h"
#include "io.h"
#include "box.h"
#include "gid.h"

// ToDo:
//
// numberRecords isn't known until we look the data (and requires a
// reduce).  Do we really need to set it?
//   Yes. And now we are setting it. / tomaso
//
// possibly add more constraints to the subset: species, rMin, rMax
//
// Replace function declarations with header includes.

typedef struct subsetWrite_parms_st
{
   int       _nFiles;
   int       _modulus;
   int       _odd;
   int       _nIdList;
   gid_type* _idList;
   gid_type  _idMin, _idMax;
   double    _xMin,  _xMax;
   double    _yMin,  _yMax;
   double    _zMin,  _zMax;
   double    _vxMin,  _vxMax;
   double    _vyMin,  _vyMax;
   double    _vzMin,  _vzMax;
   unsigned  _rPrecision;
   unsigned  _vPrecision;
   char*     _parms_info;
   char*     _filename;
   int*      _includeSpecies;
   char*     _format;
   char*     _lengthUnit;
} SUBSET_WRITE_PARMS;


static void subsetWritePio(ANALYSIS* analysis);
static void subsetWriteOvito(ANALYSIS* analysis);
static void subsetWriteBinaryCharmm(ANALYSIS* analysis);

static int rejectParticle(STATE* state, int iParticle, SUBSET_WRITE_PARMS* parms);


SUBSET_WRITE_PARMS *subsetWrite_parms(ANALYSIS *analysis)
{
   SUBSET_WRITE_PARMS *parms = ddcMalloc(sizeof(SUBSET_WRITE_PARMS)); 
	
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system;
   OBJECT* obj = (OBJECT*) analysis;
   char buf[32];
   parms->_idList = NULL; 

   object_get(obj, "nfiles",     &parms->_nFiles,     INT, 1, "0");
   object_get(obj, "modulus",    &parms->_modulus,    INT, 1, "1");
   parms->_nIdList=object_getv(obj, "idList",    &parms->_idList,     U64, IGNORE_IF_NOT_FOUND);
   object_get(obj, "idmin",      &parms->_idMin,      U64, 1, "0");
   sprintf(buf, "%"PRIu64, gid_max);
   object_get(obj, "idmax",      &parms->_idMax,      U64, 1, buf);
   object_get(obj, "odd",        &parms->_odd,        INT, 1, "0");
   object_get(obj, "rprecision", &parms->_rPrecision, INT, 1, "4");
   object_get(obj, "vprecision", &parms->_vPrecision, INT, 1, "4");
   object_get(obj, "filename",   &parms->_filename,  STRING, 1, "subset");
   object_get(obj, "format",     &parms->_format,    STRING, 1, "pio");
   object_get(obj, "lengthUnit",  &parms->_lengthUnit, LITERAL, 1, "Ang");

   qsort(parms->_idList,parms->_nIdList,sizeof(gid_type),compareGid);
   
   char** speciesNames = NULL;
   int nsIn = object_getv(obj, "species",  (void*)&speciesNames, STRING, IGNORE_IF_NOT_FOUND);

   int nSpecies=system_getNspecies(sys);
   parms->_includeSpecies = ddcMalloc(nSpecies * sizeof(int));
   for (int ii=0; ii<nSpecies; ++ii)
      parms->_includeSpecies[ii] = 1;
   
   if (nsIn > 0 )
   {
      for (int ii=0; ii<nSpecies; ++ii)
			parms->_includeSpecies[ii] = 0;
      for (int ii=0; ii<nsIn; ++ii)
      {
			SPECIES* sPtr = species_find(NULL, speciesNames[ii]);
			int iSpecies = sPtr->index;
			parms->_includeSpecies[iSpecies] = 1; 
			ddcFree(speciesNames[ii]);
      }
   }
   ddcFree(speciesNames);
   
   if (parms->_rPrecision < 4 ) parms->_rPrecision = 4;
   if (parms->_vPrecision < 4 ) parms->_vPrecision = 4;
   if (parms->_rPrecision > 4 ) parms->_rPrecision = 8;
   if (parms->_vPrecision > 4 ) parms->_vPrecision = 8;
   
   // Compute a length that is guaranteed to be larger that the
   // longest box edge to use as the default for xMax, etc.  Use
   // Compute a length that is guaranteed to be larger that the
   // longest box edge to use as the default for xMax, etc.  Use
   // -maxSize for default value of xMin, etc.
   THREE_VECTOR bbox = box_get_boundingbox(NULL);
   double maxSize = MAX(bbox.x, bbox.y);
   maxSize = MAX(bbox.z, maxSize);
   maxSize = units_convert(maxSize, "length_internal", "l");
   sprintf(buf, "%e", -maxSize);
   object_get(obj, "xmin", &parms->_xMin, WITH_UNITS, 1, buf,"l","length_internal");
   object_get(obj, "ymin", &parms->_yMin, WITH_UNITS, 1, buf,"l","length_internal");
   object_get(obj, "zmin", &parms->_zMin, WITH_UNITS, 1, buf,"l","length_internal");
   sprintf(buf, "%e", maxSize);
   object_get(obj, "xmax", &parms->_xMax, WITH_UNITS, 1, buf,"l","length_internal");
   object_get(obj, "ymax", &parms->_yMax, WITH_UNITS, 1, buf,"l","length_internal");
   object_get(obj, "zmax", &parms->_zMax, WITH_UNITS, 1, buf,"l","length_internal");
   maxSize = DBL_MAX; 
   sprintf(buf, "%e", -maxSize);
   object_get(obj, "vxmin", &parms->_vxMin, WITH_UNITS, 1, buf,"l/t","velocity_internal");
   object_get(obj, "vymin", &parms->_vyMin, WITH_UNITS, 1, buf,"l/t","velocity_internal");
   object_get(obj, "vzmin", &parms->_vzMin, WITH_UNITS, 1, buf,"l/t","velocity_internal");
   sprintf(buf, "%e", maxSize);
   object_get(obj, "vxmax", &parms->_vxMax, WITH_UNITS, 1, buf,"l/t","velocity_internal");
   object_get(obj, "vymax", &parms->_vyMax, WITH_UNITS, 1, buf,"l/t","velocity_internal");
   object_get(obj, "vzmax", &parms->_vzMax, WITH_UNITS, 1, buf,"l/t","velocity_internal");
/*
   parms->_xMin = units_convert(parms->_xMin, "l", "length_internal");
   parms->_xMax = units_convert(parms->_xMax, "l", "length_internal");
   parms->_yMin = units_convert(parms->_yMin, "l", "length_internal");
   parms->_yMax = units_convert(parms->_yMax, "l", "length_internal");
   parms->_zMin = units_convert(parms->_zMin, "l", "length_internal");
   parms->_zMax = units_convert(parms->_zMax, "l", "length_internal");
*/
   
   char string[1024];
   double lengthConvert = units_convert(1, "length_internal", "Angstrom");
   double velocityConvert = units_convert(1, "velocity_internal", "Angstrom/fs");
   snprintf(string, 1023, "idmin = %"PRIu64"; idmax = %"PRIu64"; modulus = %d; odd = %d;\n"
	    "xmin = %f Ang; xmax = %f Ang;\n"
	    "ymin = %f Ang; ymax = %f Ang;\n"
	    "zmin = %f Ang; zmax = %f Ang;\n"
	    "vxmin = %f Ang/fs; vxmax = %f Ang/fs;\n"
	    "vymin = %f Ang/fs; vymax = %f Ang/fs;\n"
	    "vzmin = %f Ang/fs; vzmax = %f Ang/fs;\n",
	    parms->_idMin, parms->_idMax, parms->_modulus, parms->_odd,
	    parms->_xMin*lengthConvert, parms->_xMax*lengthConvert,
	    parms->_yMin*lengthConvert, parms->_yMax*lengthConvert,
	    parms->_zMin*lengthConvert, parms->_zMax*lengthConvert,
	    parms->_vxMin*velocityConvert, parms->_vxMax*velocityConvert,
	    parms->_vyMin*velocityConvert, parms->_vyMax*velocityConvert,
	    parms->_vzMin*velocityConvert, parms->_vzMax*velocityConvert);
   parms->_parms_info = strdup(string); 
   return parms;
}

void subsetWrite(ANALYSIS* analysis)
{
   SUBSET_WRITE_PARMS* parms = (SUBSET_WRITE_PARMS *)analysis->parms; 
   if (! strcmp(parms->_format, "pio")) subsetWritePio(analysis);
   if (! strcmp(parms->_format, "ovito")) subsetWriteOvito(analysis);
   if (! strcmp(parms->_format, "binaryCharmm")) subsetWriteBinaryCharmm(analysis);
}

      
void subsetWritePio(ANALYSIS* analysis)
{
   SUBSET_WRITE_PARMS* parms = (SUBSET_WRITE_PARMS *)analysis->parms; 
   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state; 
   SPECIES** species = state->species;
   GROUP** group = state->group;
   gid_type* label = state->label;
   unsigned nlocal = simulate->system->nlocal;

   int nfields = 11;
   char fieldTypes[nfields*3+1]; // nfields*3 + 1
   char fieldNames[] = "checksum id pinfo rx  ry  rz  vx  vy  vz U domain";
   char fieldFormat[] = { 'u', 'b', 'b', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'b'};
   int  fieldSize[]   = {  4,   8,   1,   4,   4,   4,   4,   4,   4,   4,   8};

	// Ascii stuff {
	/*
	  PioSet(file, "field_names", "checksum id class type group rx ry rz vx vy vz");
	  PioSet(file, "field_types", "u u s s s f f f f f f" );
	  const char fmt[] = "%08x %12llu %s %s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e";
	*/
	// }
   
   PINFO_CODEC* pinfo = pinfoEncodeInit();
   unsigned pinfoFieldSize = bFieldSize(pinfoMaxIndex(pinfo));
   unsigned gidFieldSize = bFieldSize(mpiMaxVal(label, nlocal,COMM_LOCAL));
   unsigned domainFieldSize = bFieldSize(getSize(0));
   
   fieldSize[1] =  gidFieldSize;
   fieldSize[2] =  pinfoFieldSize;
   fieldSize[10] = domainFieldSize;
   if (parms->_rPrecision == 8)
      fieldSize[3] = fieldSize[4] = fieldSize[5] = 8;
   if (parms->_vPrecision == 8)
      fieldSize[6] = fieldSize[7] = fieldSize[8] = 8;
   
   for (int ii=0; ii<nfields; ++ii)
      snprintf(fieldTypes+3*ii, 4, " %.1s%1d", fieldFormat+ii, fieldSize[ii]);
   
   unsigned* offset = ddcMalloc(nfields*sizeof(unsigned));
   offset[0] = 0;
   int lrec = fieldSize[0];
   for (int ii=1; ii<nfields; ++ii)
   {
      lrec += fieldSize[ii];
      offset[ii] = offset[ii-1] + fieldSize[ii-1];
   }
   
   CreateSnapshotdir(simulate, NULL);
   char filename[1024];
   snprintf(filename, 1023,"%s/%s", simulate->snapshotdir, parms->_filename);
   PFILE* file = Popen(filename, "w", COMM_LOCAL);
   PioSet(file, "recordLength", lrec);
   PioSet(file, "datatype", FIXRECORDBINARY);


	/* Precount all records, to get a valid header. */
	gid_type nRecords = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
      if (! rejectParticle(state, ii, parms))
         ++nRecords;
	{
	  gid_type tmp = nRecords;
	  MPI_Allreduce(&tmp, &nRecords, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	}
   PioSet(file, "numberRecords", nRecords); 
   PioSet(file, "checksum", CRC32);
   PioSet(file, "nfields", nfields);
   PioSet(file, "field_names", fieldNames );
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "misc_info", parms->_parms_info);

   PioSet(file, "misc_info", "groups =");
   for (unsigned ii=0; ii<pinfo->_nGroups; ++ii) PioSet(file, "misc_info", pinfo->_gNames[ii]);
   PioSet(file, "misc_info", ";");

   PioSet(file, "misc_info", "species =");
   for (unsigned ii=0; ii<pinfo->_nSpecies; ++ii) PioSet(file, "misc_info", pinfo->_sNames[ii]);
   PioSet(file, "misc_info", ";");

   PioSet(file, "misc_info", "ntypes =");
   for (unsigned ii=0; ii<pinfo->_nTypes; ++ii) PioSet(file, "misc_info", pinfo->_tNames[ii]);
   PioSet(file, "misc_info", ";");

   if (parms->_nFiles > 0) PioSet(file, "ngroup", parms->_nFiles);

   if (getRank(0) == 0) write_fileheader(file, simulate, "subset" );

   unsigned char line[1024];
   double length_convert = units_convert(1.0,NULL,"l"); 
   double velocity_convert = units_convert(1.0,NULL,"velocity"); 
   double energy_convert = units_convert(1.0,NULL,"eV"); 
   gid_type nWritten = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (rejectParticle(state, ii, parms)) continue;

		++nWritten;

      const gid_type gid = state->label[ii];
      const double rx = state->rx[ii];
      const double ry = state->ry[ii];
      const double rz = state->rz[ii];
      const double vx = state->vx[ii];
      const double vy = state->vy[ii];
      const double vz = state->vz[ii];
    //const int iSpecies = state->species[ii]->index;
      const double U = state->potentialEnergy[ii];

      double f8;
      float f4;
      LONG64 pinfoIndex = pinfoEncode(group[ii], species[ii], pinfo);
      bFieldPack(line+offset[1], gidFieldSize,   gid);
      bFieldPack(line+offset[2], pinfoFieldSize, pinfoIndex);
      if (parms->_rPrecision == 8)
      {
         f8 = rx*length_convert;   copyBytes(line+offset[3], &f8, 8);
         f8 = ry*length_convert;   copyBytes(line+offset[4], &f8, 8);
         f8 = rz*length_convert;   copyBytes(line+offset[5], &f8, 8);
      }
      else
      {
         f4 = rx*length_convert;   copyBytes(line+offset[3], &f4, 4);
         f4 = ry*length_convert;   copyBytes(line+offset[4], &f4, 4);
         f4 = rz*length_convert;   copyBytes(line+offset[5], &f4, 4);
      }
      if (parms->_vPrecision == 8)
      {
         f8 = vx*velocity_convert;   copyBytes(line+offset[6], &f8, 8);
         f8 = vy*velocity_convert;   copyBytes(line+offset[7], &f8, 8);
         f8 = vz*velocity_convert;   copyBytes(line+offset[8], &f8, 8);
      }
      else
      {
         f4 = vx*velocity_convert;   copyBytes(line+offset[6], &f4, 4);
         f4 = vy*velocity_convert;   copyBytes(line+offset[7], &f4, 4);
         f4 = vz*velocity_convert;   copyBytes(line+offset[8], &f4, 4);
      }
      f4 = U*energy_convert;    copyBytes(line+offset[9], &f4, 4);
      bFieldPack(line+offset[10], domainFieldSize, getRank(0));
      
      int i4 =  checksum_crc32_table(line+offset[1],lrec-offset[1]);
      copyBytes(line, &i4, 4);
      Pwrite(line, lrec, 1, (PFILE *)file);
   }

   pinfoFree(pinfo);
   Pclose(file);
   ddcFree(offset);

	// Check that pre-counting yields number of records written
	{
	  gid_type tmp = nWritten;
	  MPI_Allreduce(&tmp, &nWritten, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	  if(getRank(0) == 0)
		 assert(nRecords == nWritten && "Different filters in counting and writing phase. Go fix!" != NULL);
	  MPI_Barrier(COMM_LOCAL);
	}
}

void subsetWriteOvito(ANALYSIS* analysis)
{
   SUBSET_WRITE_PARMS* parms = (SUBSET_WRITE_PARMS *)analysis->parms; 
   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state; 
   unsigned nlocal = simulate->system->nlocal;

   char fieldNames[] = "gid rx  ry  rz  vx  vy  vz species U";
   const char fmt[] = "%14"PRIu64" %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e %4d %21.13e\n";

   CreateSnapshotdir(simulate, NULL);
   char filename[1024];
   snprintf(filename, 1023,"%s/%s", simulate->snapshotdir, parms->_filename);
   PFILE* file = Popen(filename, "w", COMM_LOCAL);

	// Count the number of atoms we will write (global sum)
	gid_type nRecords = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
      if (! rejectParticle(state, ii, parms))
         nRecords++;
   {
      gid_type tmp = nRecords;
      MPI_Allreduce(&tmp, &nRecords, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	}

   PioSet(file, "ngroup", 1);

   if (getRank(0) == 0)
   {
      Pprintf(file, "%d\n", nRecords);
      Pprintf(file, "%s\n", fieldNames);
   }
   
   double length_convert = units_convert(1.0,NULL,"l"); 
   double velocity_convert = units_convert(1.0,NULL,"velocity"); 
   double energy_convert = units_convert(1.0,NULL,"eV"); 
   gid_type nWritten = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (rejectParticle(state, ii, parms)) continue;

      ++nWritten;

      gid_type gid = state->label[ii];
      double rx = state->rx[ii] * length_convert;
      double ry = state->ry[ii] * length_convert;
      double rz = state->rz[ii] * length_convert;
      double vx = state->vx[ii] * velocity_convert;
      double vy = state->vy[ii] * velocity_convert;
      double vz = state->vz[ii] * velocity_convert;
      const int iSpecies = state->species[ii]->index;
      double U = state->potentialEnergy[ii] * energy_convert;
      
      Pprintf(file, fmt, gid, rx, ry, rz, vx, vy, vz, iSpecies, U);
   }
   
   Pclose(file);
   
	// Check that pre-counting yields number of records written
	{
	  gid_type tmp = nWritten;
	  MPI_Allreduce(&tmp, &nWritten, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	  if(getRank(0) == 0)
		 assert(nRecords == nWritten && "Different filters in counting and writing phase. Go fix!" != NULL);
	  MPI_Barrier(COMM_LOCAL);
	}
}
void subsetWriteBinaryCharmm(ANALYSIS* analysis)
{
   SUBSET_WRITE_PARMS* parms = (SUBSET_WRITE_PARMS *)analysis->parms; 
   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state; 
   unsigned nlocal = simulate->system->nlocal;
   PINFO_CODEC* pinfo = pinfoEncodeInit();

   unsigned pinfoFieldSize = bFieldSize(pinfoMaxIndex(pinfo));
   assert(pinfoFieldSize <= 4); 
   int nfields = 5;
   char fieldTypes[nfields*3+1]; // nfields*3 + 1
   char fieldNames[] = "id pinfo rx ry  rz";
   char fieldFormat[] = {  'u', 'u', 'f', 'f', 'f'};
   int  fieldSize[]   = {   8,   4,   4,   4,   4 };
   int lfieldUnits = 3*(strlen(parms->_lengthUnit) +1) + 4;
   char fieldUnits[lfieldUnits];
   snprintf(fieldUnits, lfieldUnits, "1 1 %s %s %s", parms->_lengthUnit, parms->_lengthUnit, parms->_lengthUnit);
   
   for (int ii=0; ii<nfields; ++ii)
      snprintf(fieldTypes+3*ii, 4, " %.1s%1d", fieldFormat+ii, fieldSize[ii]);
   
   unsigned* offset = ddcMalloc(nfields*sizeof(unsigned));
   int lrec = fieldSize[0];
   offset[0] = 0;
   for (int ii=1; ii<nfields; ++ii)
   {
      lrec += fieldSize[ii];
      offset[ii] = offset[ii-1] + fieldSize[ii-1];
   }
   
   CreateSnapshotdir(simulate, NULL);  //   CHECK
   char filename[1024];
   snprintf(filename, 1023,"%s/%s", simulate->snapshotdir, parms->_filename); //CHECK
   PFILE* file = Popen(filename, "w", COMM_LOCAL);

	/* Precount all records, to get a valid header. */
	gid_type nRecords = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
      if (! rejectParticle(state, ii, parms)) ++nRecords;
	{
	  gid_type tmp = nRecords;
	  MPI_Allreduce(&tmp, &nRecords, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	}
   PioSet(file, "recordLength", lrec);
   PioSet(file, "datatype", FIXRECORDBINARY);
   PioSet(file, "numberRecords", nRecords); 
   PioSet(file, "checksum", NONE);
   PioSet(file, "nfields", nfields);
   PioSet(file, "field_names", fieldNames );
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "field_units", fieldUnits);
  if (parms->_nFiles > 0) PioSet(file, "ngroup", parms->_nFiles);


	PioSet(file, "misc_info", "random = NONE;\n");

	PioSet(file, "misc_info", "nrandomFieldSize = 0;\n");

   PioSet(file, "misc_info", "types =");
   for (unsigned ii=0; ii<pinfo->_nTypes; ++ii) PioSet(file, "misc_info", pinfo->_tNames[ii]);
   PioSet(file, "misc_info", ";\n");

   PioSet(file, "misc_info", "groups =");
   for (unsigned ii=0; ii<pinfo->_nGroups; ++ii) PioSet(file, "misc_info", pinfo->_gNames[ii]);
   PioSet(file, "misc_info", ";\n");

   PioSet(file, "misc_info", "species =");
   for (unsigned ii=0; ii<pinfo->_nSpecies; ++ii) PioSet(file, "misc_info", pinfo->_sNames[ii]);
   PioSet(file, "misc_info", ";\n");

   PioSet(file, "misc_info", parms->_parms_info);


   if (getRank(0) == 0) write_fileheader(file, simulate, "subset" );

   unsigned char line[1024];
   double cL = units_convert(1.0,NULL,parms->_lengthUnit); 
   gid_type nWritten = 0;
   THREE_VECTOR corner = box_get_corner(NULL); 
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (rejectParticle(state, ii, parms)) continue;

		++nWritten;

      const gid_type gid = state->label[ii];
      const double rx = state->rx[ii]-corner.x;
      const double ry = state->ry[ii]-corner.y;
      const double rz = state->rz[ii]-corner.z;
      const unsigned index = (unsigned)pinfoEncode(state->group[ii], state->species[ii], pinfo);

      float f4;
                       copyBytes(line+offset[0],&gid, 8);
                       copyBytes(line+offset[1],&index, 4);
         f4 = rx*cL;   copyBytes(line+offset[2], &f4, 4);
         f4 = ry*cL;   copyBytes(line+offset[3], &f4, 4);
         f4 = rz*cL;   copyBytes(line+offset[4], &f4, 4);
      Pwrite(line, lrec, 1, (PFILE *)file);
   }
   pinfoFree(pinfo);
   Pclose(file);
   ddcFree(offset);

	// Check that pre-counting yields number of records written
	{
	  gid_type tmp = nWritten;
	  MPI_Allreduce(&tmp, &nWritten, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
	  if(getRank(0) == 0)
		 assert(nRecords == nWritten && "Different filters in counting and writing phase. Go fix!" != NULL);
	  MPI_Barrier(COMM_LOCAL);
	}
}


void subsetWrite_eval(ANALYSIS* analysis)
{
   // empty
}

// returns non-zero if particle is rejected.
// return zero if particle is in the subset.
int rejectParticle(STATE* state, int iParticle, SUBSET_WRITE_PARMS* parms)
{
   const gid_type gid = state->label[iParticle];
   const double rx = state->rx[iParticle];
   const double ry = state->ry[iParticle];
   const double rz = state->rz[iParticle];
   const double vx = state->vx[iParticle];
   const double vy = state->vy[iParticle];
   const double vz = state->vz[iParticle];
   const int iSpecies = state->species[iParticle]->index;

   // conditionals to *reject* partices
   if (gid < parms->_idMin )        return 1;
   if (gid > parms->_idMax )        return 1;
   if (gid % parms->_modulus != 0 ) return 1;
   if (parms->_odd && gid % 2 ==0 ) return 1;
   if (rx > parms->_xMax )          return 1;
   if (ry > parms->_yMax )          return 1;
   if (rz > parms->_zMax )          return 1;
   if (rx < parms->_xMin )          return 1;
   if (ry < parms->_yMin )          return 1;
   if (rz < parms->_zMin )          return 1;
   if (vx > parms->_vxMax )          return 1;
   if (vy > parms->_vyMax )          return 1;
   if (vz > parms->_vzMax )          return 1;
   if (vx < parms->_vxMin )          return 1;
   if (vy < parms->_vyMin )          return 1;
   if (vz < parms->_vzMin )          return 1;
   if (parms->_includeSpecies[iSpecies] == 0) return 1;
   if (parms->_idList != NULL && bsearch(&gid,parms->_idList,parms->_nIdList,sizeof(gid_type),compareGid) == NULL) return 1;

   return 0;
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
