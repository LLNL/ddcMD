#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define  _LARGE_FILE
#include "collection.h"

#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "object.h"
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "pioVariableRecordHelper.h"
#include "pinfo.h"
#include "ioUtils.h"
#include "utilities.h"
#include "ddcMalloc.h"
#include "crc32.h"
#include "units.h"
#include "mpiUtils.h"
#include "system.h"
#include "state.h"
#include "mpiUtils.h"
#include "check_line.h"
#include "external.h"



extern FILE* ddcfile;

static void validateFieldMap(enum IO_FIELDS* fieldMap, unsigned nFields);
static unsigned findIof(int target, const enum IO_FIELDS* array, unsigned size);


#define UNASSIGNED 0xffffffffffffffff


void collection_readVARRECORDASCII(COLLECTION* c, PFILE* pfile)
{
   readline_init(pfile);
   PIO_VARIABLE_RECORD_ASCII_HELPER* helper = (PIO_VARIABLE_RECORD_ASCII_HELPER*) pfile->helper;
   unsigned nRecords = pvrah_nRecords(pfile->buf, pfile->bufsize, helper->delimiter);
	resize(c->size + nRecords, 2, c->state);
   collection_readASCII(c->state, nRecords, pfile);
   c->size += nRecords;
   timestampBarrier("Finished collection_readVARRECORDASCII\n", COMM_LOCAL);
}

void collection_readFIXRECORDASCII(COLLECTION* c, PFILE* pfile)
{
   readline_init(pfile);
   PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) pfile->helper;
   unsigned nRecords = pfile->bufsize/helper->lrec;
   assert(pfile->bufsize%helper->lrec == 0);
	resize(c->size + nRecords, 2, c->state);
   collection_readASCII(c->state, nRecords, pfile);
   c->size += nRecords;
 
   timestampBarrier("Finished collection_readFIXRECORDASCII\n", COMM_LOCAL);
}
void collection_readFIXRECORDBINARY(COLLECTION* c, PFILE* pfile)
{
	PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) pfile->helper;
	unsigned recordLength = helper->lrec;
	unsigned nRecords = pfile->bufsize/recordLength;
	assert(pfile->bufsize%recordLength == 0);
	resize(c->size + nRecords, 2, c->state);
   collection_readBINARY(c->state, recordLength, nRecords, pfile);
   c->size += nRecords;
}

enum hasField {No,Yes,Unknown};
/*
 Example:
  *  char* tmp[3];
  *   *  int n = object_get(obj, "foo", &tmp, STRING, 3, "NONE");
  *    *  for (int ii=0; ii<n; ++ii)
  *     *  {
  *      *      foo_init(tmp[ii]);
  *       *      ddcFree(tmp[ii]);
  *        *   }
  *        
  *        */
void collection_readASCII(STATE *state, int size, PFILE *pfile)
{
	char  message[256], *class, *tail;
	sprintf(message, "Processor %d: Is Starting to read %d lines of ASCII format file", getRank(0), size);
	timestamp(message);
	long msgInterval = lrint(pow(10.0,round(log(size/10.0)/log(10.0))));
	if (msgInterval <=  0 ) msgInterval = size; 
	
	double length_convert = units_convert(1.0,"l",NULL); 
	double time_convert = units_convert(1.0,"t",NULL); 
	double velocity_convert = length_convert/time_convert; 
   int maxRecordLength=2048; 
   char record[maxRecordLength];
	int checksum_ok;
	int error=0; 	
   OBJECT* hobj = pfile->headerObject;
   char *string; 
   object_get(hobj, "random", &string, STRING, 1, "NotSet");  //NotSet can be used as a test to check in file use  old format for random number portion of the input line
   int hasRandomField = Yes; 
   if (strcmp(string,"NotSet")==0) hasRandomField=Unknown;
   if (strcmp(string,"NONE")==0)   hasRandomField=No;
   int useRandomDefault = 0; 
   RANDOM *random = system_getRandom(NULL); 
   if (random != NULL && hasRandomField == No)  useRandomDefault = 1; 

   char labelFmt[16];  
   strncpy(labelFmt,"%"PRIu64" %ln",16);  //default format; 
   if (object_testforkeyword(hobj, "field_format"))
   {
      char *format[2]; 
      int nFormat=object_get(hobj, "field_format", &format, STRING, 2, "NotSet"); 
      if (format[1] != NULL) 
      {
         char end = format[1][strlen(format[1])-1];
         assert( end == 'x' || end =='u');
         if (end == 'x') strncpy(labelFmt,"%"PRIx64" %ln",15); 
         ddcFree(format[1]); 
      }
      ddcFree(format[0]); 
   }
	for (int ii = 0; ii < size; ii++)
	{
	   gid_type label; 
		//char atomname[16],groupname[16],sep[16];
      Pfgets(record, maxRecordLength-1, pfile);
      size_t recordLength = strnlen(record,maxRecordLength);
      pfile->recordIndex++; 
      int rc; 
      char *line = checkRecord(&rc,recordLength,record,pfile); 
		checksum_ok = !(rc & 2) ;
		checksum_ok = 1 ;
		parsehead(line,labelFmt,&label,&class,&tail,checksum_ok); 
      char *end; 
      char* atomname = strtok_r(tail, " ",&end);
      SPECIES *species = species_find(NULL,atomname); 
      char* groupname= strtok_r(NULL, " ",&end);
      assert(groupname+strlen(groupname)+1==end); 
      GROUP *group = group_find(NULL,groupname); 
      double rx       = strtod(end, &end);
      double ry       = strtod(end, &end);
      double rz       = strtod(end, &end);
      double vx       = strtod(end, &end);
      double vy       = strtod(end, &end);
      double vz       = strtod(end, &end);

		if (state->label != NULL) state->label[ii] = label; 
		if (state->species != NULL) state->species[ii] = species;
		if (state->group != NULL) state->group[ii] = group;
		if (state->rx != NULL) state->rx[ii] =   length_convert*rx;
		if (state->ry != NULL) state->ry[ii] =   length_convert*ry;
		if (state->rz != NULL) state->rz[ii] =   length_convert*rz;
		if (state->vx != NULL) state->vx[ii] = velocity_convert*vx;
		if (state->vy != NULL) state->vy[ii] = velocity_convert*vy;
		if (state->vz != NULL) state->vz[ii] = velocity_convert*vz;

      if (hasRandomField != No && random != NULL) 
      {
           void *randomParms = random_getParms(random, ii);
           char *start = end; 
           end = random->parse(start, randomParms);
           if (end == start ) useRandomDefault=1; 
      }
		if (state->group != NULL) 
      {
         if (group->parse != NULL) 
         {
            char *start=index(end,'|'); // if no "|" then end == start
            if (start != NULL) 
            {
               start+=1; 
               end = group->parse(group, start, ii); // if parse fails end == start
            }
            if (end == start) group->defaultValue(group, label, ii);   //if end == start there is not a good field to parse. use default value. 
         }
      }
      if ( (state->species != NULL && species == NULL) || (state->species!=NULL && group == NULL))
      {
         printf( "Problem detected reading data file on task %d\n"
               "   atomname = %s,  groupname = %s\n   tail = X%sX\n   ii = %d\n", 
               getRank(0), atomname, groupname, tail, ii);
      }
      if (state->atomtype != NULL) state->atomtype[ii] = state->species[ii]->index + (state->group[ii]->index << 16);
      if ((ii + 1)%msgInterval == 0)
      {
         sprintf(message, "Processor %d: Finished Reading line %d of ASCII format file", getRank(0), ii + 1);
         timestamp(message);
      }
   }
   if (random != NULL) random->useDefault = useRandomDefault; 
   sprintf(message, "Processor %d: Finished Reading line %d of ASCII format file", getRank(0), size);
   timestamp(message);
   timestamp("Finished Reading ASCII format file");
   WAIT(0);
   if ( error ) { abortAll(1); } 
}
void collection_readBINARY(STATE *state, int recordLength, int nRecords,  PFILE* pfile)
{
   OBJECT* hobj = pfile->headerObject;
   char* string;
   object_get(hobj, "datatype", &string, STRING, 1, " ");
   assert(strcasecmp(string, "FIXRECORDBINARY") == 0);
   enum PIO_ENUMS checksum_type=PIO_NONE;
   object_get(hobj, "checksum", &string, STRING, 1, "NONE");
   if (strcasecmp(string, "CRC32") == 0) checksum_type=CRC32;

   RANDOM *random = system_getRandom(NULL); 
   int randomFieldSize=0; 
   object_get(hobj, "random", &string, STRING, 1, "NotSet");  //NotSet can be used as a test to check in file use  old format for random number portion of the input line
   object_get(hobj, "randomFieldSize",&randomFieldSize,INT,1,"0"); 
   int hasRandomField=No;
   if (strcmp(string,"NotSet")==0  ) 
   {
      hasRandomField=No;
      randomFieldSize =0; 
      if (random!=NULL) 
      {
         //hasRandomField=Yes;
         hasRandomField=No;
        void *randomParms = random_getParms(random, 0);
        randomFieldSize = random->bwrite(NULL,randomParms); 
      }
   }
   else
   {
      hasRandomField=Yes;
      if (strcmp(string,"NONE")==0)   
      {
         hasRandomField=No;
         randomFieldSize=0; 
      }
   }
   if (random != NULL && hasRandomField == No) random->useDefault = 1; 
   ddcFree(string);
   unsigned checksum=0;
   unsigned key;
   object_get(hobj, "endian_key", &key, INT, 1, "0");
   ioUtils_setSwap(key);

   PINFO_CODEC* pinfo = pinfoDecodeInit(hobj);
   assert(pinfo != NULL);
   char** fieldTypes = NULL;
   unsigned nFields = object_getv(hobj, "field_types", (void*)&fieldTypes, STRING,ABORT_IF_NOT_FOUND);
   assert(nFields == pfile->nfields);
   char** fieldNames = NULL;
   nFields = object_getv(hobj, "field_names", (void*)&fieldNames, STRING,ABORT_IF_NOT_FOUND);
   assert(nFields == pfile->nfields);

   unsigned* offsetField = makeOffsets(fieldTypes, nFields);
   enum IO_FIELDS* fieldMap = makeFieldMap(fieldNames, nFields);
   validateFieldMap(fieldMap, nFields);
   unsigned char* buf = ddcMalloc(recordLength*sizeof(char));

   // initialize data that might be absent in the input file
   for (int ii=0; ii<nRecords; ++ii)
   {
      if (state->group != NULL) state->group[ii] = group_by_index(NULL, 0);
      if (state->species != NULL) state->species[ii] = species_by_index(NULL, 0);
      if (state->atomtype != NULL) state->atomtype[ii] = 0;
      if (state->vx != NULL) state->vx[ii] = 0;
      if (state->vy != NULL) state->vy[ii] = 0;
      if (state->vz != NULL) state->vz[ii] = 0;
   }
   double length_convert = units_convert(1.0,"l",NULL); 
   double time_convert = units_convert(1.0,"t",NULL); 
   double velocity_convert = length_convert/time_convert; 
   for (int ii=0; ii<nRecords; ++ii)
   {
      Pread(buf, recordLength, 1, pfile);
      if (checksum_type == CRC32) checksum = checksum_crc32(buf+4, recordLength-4);
       char* type;
      SPECIES *species=NULL; 
      GROUP *group=NULL; 
      for (unsigned jj=0; jj<nFields; ++jj)
      {
         unsigned long long tmp =  UNASSIGNED;
         unsigned char* fPtr = buf+offsetField[jj];
         switch (fieldMap[jj])
         {
            case IOF_CHECKSUM:
               tmp = mkInt(fPtr, fieldTypes[jj]);
               if (tmp != checksum)
               {
                  printf("ERROR: Checksum Error on task %d\n", getRank(0));
                  abortAll(45);
               }
               break;
            case IOF_ID:
               if (state->label != NULL) state->label[ii] = mkInt(fPtr, fieldTypes[jj]);
               break;
            case IOF_SPECIES:
               tmp = mkInt(fPtr, fieldTypes[jj])*pinfo->_nGroups;
            case IOF_PINFO:
               if (tmp == UNASSIGNED) tmp = mkInt(fPtr, fieldTypes[jj]);
               int rc = pinfoDecode(tmp, &group, &species, &type, pinfo);
               assert(strcmp(type, "ATOM") == 0);
               if (rc != 0)
               {
                  printf("ERROR: Problem decoding species/group on task %d.\n" "       data=%lld\n" "       ii=%d\n", getRank(0), tmp, ii);
                  abortAll(45);
               }
               if (state->atomtype != NULL) state->atomtype[ii] = species->index + (group->index << 16);
               if (state->species != NULL) state->species[ii] = species; 
               if (state->group   != NULL) state->group[ii] = group; 
               break;
            case IOF_RX:
               if (state->rx != NULL) state->rx[ii] = length_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_RY:
               if (state->ry != NULL) state->ry[ii] = length_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_RZ:
               if (state->rz != NULL) state->rz[ii] = length_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_VX:
               if (state->vx != NULL) state->vx[ii] = velocity_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_VY:
               if (state->vy != NULL) state->vy[ii] = velocity_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_VZ:
               if (state->vz != NULL) state->vz[ii] = velocity_convert*mkDouble(fPtr, fieldTypes[jj]);
               break;
            case IOF_NONE:
               break;
            default:
               assert(1==0); // this can't happen
         }
      }
      int offset = offsetField[nFields]; 
      if (hasRandomField == Yes) 
      {
         if (random != NULL) 
         {
            void *randomParms = random_getParms(random, ii); 
            random->bread(buf+offset, randomParms);
         }
         offset += randomFieldSize; 
      }
      if (group != NULL) 
      {
         if (group->bread != NULL) group->bread(buf+offset, group, ii);
      }
   }

   ddcFree(buf);
   ddcFree(fieldMap);
   ddcFree(offsetField);
   for (unsigned ii=0; ii<nFields; ++ii)
   {
      ddcFree(fieldTypes[ii]);
      ddcFree(fieldNames[ii]);
   }
   ddcFree(fieldNames);
   ddcFree(fieldTypes);
   pinfoFree(pinfo);
}

unsigned findIof(int target, const enum IO_FIELDS* array, unsigned size)
{
   for (unsigned ii=0; ii<size; ++ii) if ((int)array[ii] == target) return ii;
   return size;
}

/** This routine calls MPI_ABORT when the fieldMap does not contain all
 * of the necessary fields to initialize a simulation. */
void validateFieldMap(enum IO_FIELDS* fieldMap, unsigned nFields)
{
   return; //DEBUG 
   if (getRank(0) != 0)
      return;

   int error = 0;

   // The global id and the atom positions must be present.
   if (findIof(IOF_ID, fieldMap, nFields) == nFields)
      error |= 1;
   if (findIof(IOF_RX, fieldMap, nFields) == nFields)
      error |= 2;
   if (findIof(IOF_RY, fieldMap, nFields) == nFields)
      error |= 4;
   if (findIof(IOF_RZ, fieldMap, nFields) == nFields)
      error |= 8;

   // Print warning if the velocities aren't available.
   if (findIof(IOF_VX, fieldMap, nFields) == nFields ||
         findIof(IOF_VY, fieldMap, nFields) == nFields ||
         findIof(IOF_VZ, fieldMap, nFields) == nFields )
      printf("\nWARNING:  One or more velocity components could not\n"
            "be found when reading collection data from binary file.\n"
            "The missing components will be set to zero.\n\n");


   // Species and group information can be absent only if there is only
   // one possibility and can therefore be synthesized.
   if (findIof(IOF_PINFO, fieldMap, nFields) == nFields)
   {
      int nGroups;
      int nSpecies;
      group_get(NULL, NGROUPS, (void*) &nGroups);
      nSpecies= system_getNspecies(NULL);
      if (nGroups == 1 && nSpecies == 1)
         printf("\nWARNING:  Group and species information was not\n"
               "found when reading binary file.  Since only one\n"
               "group and species are specified for this simulation\n"
               "all atoms will be set to that species and group.\n\n");
      else
         error |= 16;
   }


   if (error != 0)
   {
      printf("\nERROR while reading collection from binary file.\n"
            "The following required fields are missing:\n");
      if (error &  1)  printf("   global id\n");
      if (error &  2)  printf("   rx\n");
      if (error &  4)  printf("   ry\n");
      if (error &  8)  printf("   rz\n");
      if (error & 16)  printf("   group/species info\n");
      printf("\n");

      abortAll( -1);
   }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
