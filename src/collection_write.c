#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define  _LARGE_FILE

#include "collection_write.h"

#include <string.h>
#include <mpi.h>
#include <assert.h>

#include "pio.h"
#include "simulate.h"
#include "error.h"
#include "crc32.h"
#include "ioUtils.h"
#include "pinfo.h"
#include "ddcMalloc.h"
#include "utilities.h"
#include "orderSH.h"
#include "units.h"
#include "mpiAlgorithm.h"
#include "io.h"
#include "preduce.h"
#include "format.h"

int getRank(int);

#define MAXLREC   1024


static int writeMode = WHEADER;
static int _nCopies=1;

void collection_writeMode(int mode)
{
   writeMode = mode;
   if (mode == WHEADER)
      _nCopies = 1;
}

void collection_nCopies(int n)
{
   _nCopies = n;
}
char *randomWrite(RANDOM *random, int i)
{
   void *randomParms = random_getParms(random, i); 
   return random->write(randomParms);
}
unsigned randomBwrite(unsigned char* buf, RANDOM *random, int i)
{
   void *randomParms = random_getParms(random, i); 
   unsigned size = random->bwrite(buf, randomParms);
   return size; 

}
void collection_writeBLOCK(SIMULATE*simulate, PFILE*file)
{
	STATE *state;
	SPECIES **species;
	GROUP **group;
	unsigned checksum;
	gid_type *label;
	static gid_type nglobal;
	double *rx, *ry, *rz, *vx, *vy, *vz;
	char *string, tmp_string[16],line[MAXLREC], errmsg[128]; 
   //char fmt[] = "%08x %12llu %s %s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e";
   char fmt[128]; 
   sprintf(fmt,"%s %s %s","%08x ",gidFormat()," %s %s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e");
   SYSTEM *sys = simulate->system; 
	COLLECTION *c = sys->collection;
	state = c->state;
	rx = state->rx;	ry = state->ry;	rz = state->rz;
	vx = state->vx;	vy = state->vy;	vz = state->vz;
	label = state->label;
	species = state->species;
	group = state->group;
	c->size = state->nlocal;
	nglobal = simulate->system->nglobal;
	unsigned nlocal = simulate->system->nlocal;
   RANDOM *random = system_getRandom(sys); 
   int length;
   int randomFieldSize=0;
//record length needed without random and group writes assume species type, species type , species name, and group name only need 1 byte
	int lrec = snprintf(line, MAXLREC, fmt,0,0," "," "," ",0.0,0.0,0.0,0.0,0.0,0.0 ); 
	species_get(NULL, SPECIESMAXTYPELENGTH, (void *)&length); lrec += (length-1);    // Correct for species type length larger than 1
	species_get(NULL, SPECIESMAXNAMELENGTH, (void *)&length); lrec += (length-1);   // Correct for species name length larger than 1
	group_get(NULL, GROUPMAXNAMELENGTH, (void *)&length);   ; lrec += (length-1);    // Correct for group name  length larger than 1
   if (sys->random != NULL) randomFieldSize=strlen(random->write(NULL));          lrec += randomFieldSize + 1; 
	group_get(NULL, GROUPMAXWRITELENGTH, (void *)&length);   if (length > 0)       lrec += length + 3;  // need remove for delimiter. 
   lrec++;   // increment for terminator. 
	lrec = 8*((lrec+7)/8); // pad to multiples of 8 bytes. 
	if (lrec > MAXLREC)
	{
		snprintf(errmsg, 128, "Record Length=%d exceeds MAXLREC=%d program will abort", lrec, MAXLREC);
		error_action(errmsg, ERROR_IN("collection_writeBLOCK", ABORT));
	}

	if (writeMode == WHEADER8FOLD)  nglobal *= 8;
	if (writeMode == WREPLICATE) nglobal *= _nCopies;
	PioSet(file, "recordLength", lrec);
	PioSet(file, "datatype", FIXRECORDASCII);
	PioSet(file, "numberRecords", nglobal);
	PioSet(file, "checksum", CRC32);
	int nfields = 11; 
	PioSet(file, "nfields", nfields);
	PioSet(file, "field_names", "checksum id class type group rx ry rz vx vy vz");
	PioSet(file, "field_types", "u u s s s f f f f f f" );
	PioSet(file, "field_units", "1 1 1 1 1 Ang Ang Ang Ang/fs Ang/fs Ang/fs" );
	PioSet(file, "field_format", fmt );
   char *name = "NONE";
   if (random != NULL) name = random->name; 
	PioSet(file, "misc_info", "random =");
	PioSet(file, "misc_info", name);
	PioSet(file, "misc_info", ";\nrandomFieldSize =");
   char sizeString[8]; 
   sprintf(sizeString,"%d",randomFieldSize); 
	PioSet(file, "misc_info", sizeString);
	{  // limit scope of pinfo.
	   
	   // I'm creating a PINFO_CODEC here as a cheap (and dirty) way
	   // to get quick access to the the species and group lists.
	   // I could get at them directly, but this is easier and I get
	   // to copy code I already wrote for binary restarts.
	   // I want to put the lists in the header even for ascii files
	   // because it makes life easier for the visit reader.

	   PioSet(file, "misc_info", ";\ngroups =");
	   PINFO_CODEC* pinfo = pinfoEncodeInit();
	   for (unsigned ii=0; ii<pinfo->_nGroups; ++ii)
	      PioSet(file, "misc_info", pinfo->_gNames[ii]);

	   PioSet(file, "misc_info", ";\nspecies =");
	   for (unsigned ii=0; ii<pinfo->_nSpecies; ++ii)
	      PioSet(file, "misc_info", pinfo->_sNames[ii]);

	   PioSet(file, "misc_info", ";\ntypes =");
	   for (unsigned ii=0; ii<pinfo->_nTypes; ++ii)
	      PioSet(file, "misc_info", pinfo->_tNames[ii]);

	   PioSet(file, "misc_info", ";");
	   pinfoFree(pinfo);
	}
	PioReserve(file, lrec*nlocal*_nCopies + 4096);
	if (getRank(0) == 0 && writeMode != NOWHEADER) write_fileheader(file, simulate, "particle" );
	double cLen = units_convert(1.0,NULL,"l"); 
	double time_convert = units_convert(1.0,NULL,"t"); 
	double cVel = cLen/time_convert; 
	for (int i = 0; i < c->size; i++)
	{
      double x = rx[i] ; double y = ry[i]; double z = rz[i];
      backInBox(&x, &y, &z);
		int length = snprintf(line, MAXLREC, fmt, 
                      0,label[i], species[i]->type, species[i]->name, group[i]->name, 
                      x*cLen, y*cLen, z*cLen, vx[i]*cVel, vy[i]*cVel, vz[i]*cVel);
      int length0 = strlen(line); 
		if (random != NULL)
		{
			strncat(line, " ", MAXLREC - strlen(line));
			string = randomWrite(random, i);
			strncat(line, string, MAXLREC - strlen(line));
         length +=  strlen(string) +1; 
		}
      int length1 =strlen(line); 
		if (state->group[i]->write != NULL)
		{
			strncat(line, " | ", MAXLREC - strlen(line));
			string = state->group[i]->write(state->group[i], i);
			strncat(line, string, MAXLREC - strlen(line));
         length +=  strlen(string)+3; 
		}
		if (length > lrec-1)
		{
			snprintf(errmsg, 128, "For particle with label %"PRIu64" write exceeds lrec=%d lengths=%d %d %d program will abort", label[i], lrec, length0, length1,length);
			error_action(errmsg, ERROR_IN("collection_writeBLOCK", ABORT));
		}
		for (int l=length; l < lrec ; l++) line[l] = (char)' ';
		line[lrec-1] = (char)'\n';
  		checksum = checksum_crc32_table((unsigned char *)line+8,lrec-8);
		sprintf(tmp_string,"%08x",checksum);
		memcpy(line,tmp_string,8);
		Pwrite(line, lrec, 1, (PFILE *)file);
	}
	/*fflush(file);  */
}

//////////////////////////////////////////////////////////////////////
void collection_writeBLOCK_binary(SIMULATE* simulate, PFILE* file)
{
	char  errmsg[128];
	unsigned char line[MAXLREC];
   SYSTEM *sys = simulate->system; 
	COLLECTION* c = sys->collection;
	STATE* state = c->state;
	double* rx = state->rx; double* vx = state->vx;	
	double* ry = state->ry; double* vy = state->vy;	
	double* rz = state->rz; double* vz = state->vz;	
	gid_type* label = state->label;
	SPECIES** species = state->species;
	GROUP** group = state->group;
	gid_type nglobal = simulate->system->nglobal;
	unsigned nlocal = simulate->system->nlocal;
   RANDOM *random = system_getRandom(sys); 
 
	int nfields = 9; 
	char fieldTypes[28]; // nfields*3 + 1
	char fieldNames[] = "checksum id pinfo rx  ry  rz  vx  vy  vz";
	char fieldFormat[] = { 'u', 'b', 'b', 'f', 'f', 'f', 'f', 'f', 'f'};
	unsigned  fieldSize[]   = {  4,   8,   1,   8,   8,   8,   8,   8,   8};

	PINFO_CODEC* pinfo = pinfoEncodeInit();
	unsigned pinfoFieldSize = bFieldSize(pinfoMaxIndex(pinfo));
	unsigned gidFieldSize = bFieldSize(mpiMaxVal(label, nlocal,COMM_LOCAL));
	
	fieldSize[1] = gidFieldSize;
	fieldSize[2] = pinfoFieldSize;
	if (simulate->checkpointprecision == BRIEF_CHECKPOINT)
		fieldSize[6] = fieldSize[7] = fieldSize[8] = 4;
	
	for (int ii=0; ii<nfields; ++ii)
		snprintf(fieldTypes+3*ii, 4, " %.1s%1d", fieldFormat+ii, fieldSize[ii]);
	
// May need two extra offsets. offset[nfields] of the offset to write random stuff
// offset[nfields+1] is the offset to write group stuff. 
	unsigned* offset = ddcMalloc((nfields+2)*sizeof(unsigned));
	offset[0] = 0;
	unsigned lrec = fieldSize[0];
	for (int ii=1; ii<nfields; ++ii)
	{
		lrec += fieldSize[ii];
		offset[ii] = offset[ii-1] + fieldSize[ii-1];
	}
	offset[nfields] = lrec;
   unsigned randomFieldSize = 0; 
   if (random !=NULL)  
   {
      void *randomParms = random_getParms(random, 0);
      randomFieldSize = random->bwrite(NULL,randomParms); 
   }
   lrec += randomFieldSize; 
	offset[nfields+1] = lrec;
	int max;
	group_get(NULL, GROUPMAXBINARYWRITELENGTH, (void *)&max);
	lrec += max;

	if (writeMode == WHEADER8FOLD)
		nglobal *= 8;
	if (writeMode == WREPLICATE)
		nglobal *= _nCopies;
	PioSet(file, "recordLength", lrec);
	PioSet(file, "datatype", FIXRECORDBINARY);
	PioSet(file, "numberRecords", nglobal);
	PioSet(file, "checksum", CRC32);
	PioSet(file, "nfields", nfields);
	PioSet(file, "field_names", fieldNames);
	PioSet(file, "field_types", fieldTypes);
   char *name = "NONE";
   if (random != NULL) name = random->name; 
	PioSet(file, "misc_info", "random =");
	PioSet(file, "misc_info", name);
	PioSet(file, "misc_info", ";\nrandomFieldSize =");
   char sizeString[8]; 
   sprintf(sizeString,"%d",randomFieldSize); 
	PioSet(file, "misc_info", sizeString);
	PioSet(file, "misc_info", ";\ngroups =");
	for (unsigned ii=0; ii<pinfo->_nGroups; ++ii)
		PioSet(file, "misc_info", pinfo->_gNames[ii]);
	PioSet(file, "misc_info", ";\nspecies =");
	for (unsigned ii=0; ii<pinfo->_nSpecies; ++ii)
		PioSet(file, "misc_info", pinfo->_sNames[ii]);
	PioSet(file, "misc_info", ";\ntypes =");
	for (unsigned ii=0; ii<pinfo->_nTypes; ++ii)
		PioSet(file, "misc_info", pinfo->_tNames[ii]);
	PioSet(file, "misc_info", ";");
	PioReserve(file, lrec*nlocal*_nCopies + 4096);
	
	if (getRank(0) == 0 && writeMode != NOWHEADER) write_fileheader(file, simulate, "particle" );
	if (lrec > MAXLREC)
	{
		snprintf(errmsg, 128, "Record Length=%d exceeds MAXLREC=%d program will abort", lrec, MAXLREC);
		error_action(errmsg, ERROR_IN("collection_writeBLOCK_binary", ABORT));
	}
	c->size = state->nlocal;
	double cLen = units_convert(1.0,NULL,"l"); 
	double time_convert = units_convert(1.0,NULL,"t"); 
	double cVel = cLen/time_convert; 
	for (int ii=0; ii<c->size; ++ii)
	{
      for (unsigned i=0;i<lrec;i++) line[i]=0; 
		double f8;
		float f4;
		LONG64 pinfoIndex = pinfoEncode(group[ii], species[ii], pinfo);
		bFieldPack(line+offset[1], gidFieldSize,   label[ii]);
		bFieldPack(line+offset[2], pinfoFieldSize, pinfoIndex);
	   double x = rx[ii] ; double y = ry[ii]; double z = rz[ii];
	   backInBox(&x, &y, &z);
		f8 = x*cLen;   copyBytes(line+offset[3], &f8, 8);
		f8 = y*cLen;   copyBytes(line+offset[4], &f8, 8);
		f8 = z*cLen;   copyBytes(line+offset[5], &f8, 8);
		if (simulate->checkpointprecision != BRIEF_CHECKPOINT)
		{
			f8 = vx[ii]*cVel;   copyBytes(line+offset[6], &f8, 8);
			f8 = vy[ii]*cVel;   copyBytes(line+offset[7], &f8, 8);
			f8 = vz[ii]*cVel;   copyBytes(line+offset[8], &f8, 8);
		}
		else
		{
			f4 = vx[ii]*cVel;   copyBytes(line+offset[6], &f4, 4);
			f4 = vy[ii]*cVel;   copyBytes(line+offset[7], &f4, 4);
			f4 = vz[ii]*cVel;   copyBytes(line+offset[8], &f4, 4);
		}
		
		unsigned rsize = 0;
		if (random != NULL) rsize = randomBwrite(line+offset[nfields], random, ii);
      assert(rsize == randomFieldSize); 
		unsigned gsize = 0;
		if (state->group[ii]->bwrite != NULL)    
	      gsize = state->group[ii]->bwrite(line+offset[nfields+1], state->group[ii], ii);  //JNG CHECK

		if (gsize + offset[nfields+1] > lrec)
		{
			snprintf(errmsg, 128, "For particle with label %"PRIu64" lrec exceeds MAXLREC (=%d) program will abort",
				 label[ii], gsize+offset[nfields+1]);
			error_action(errmsg, ERROR_IN("collection_writeBLOCK_binary", CONTINUE));
		}
		for (unsigned jj=offset[nfields+1]+gsize; jj<lrec; ++jj ) line[jj] = '\0';

		int i4 =  checksum_crc32_table(line+offset[1],lrec-offset[1]);
		copyBytes(line, &i4, 4);
		Pwrite(line, lrec, 1, (PFILE *)file);
	}
	ddcFree(offset);
	pinfoFree(pinfo);
}

//////////////////////////////////////////////////////////////////////

void collection_writeBXYZ(SIMULATE*simulate, PFILE*file)
{
	char field_names[1024],field_types[1024];
	int i, k, n, nfields, lrec;
	unsigned mode=1; 
	double* Q = orderGetQ();
	int* C = orderGetC();
	int* Lv = orderGetLv();
	double** qnorm = orderGetqnorm();
	int nL = orderGetnL();
	STATE* state = simulate->system->collection->state;
	double* rx = state->rx;
	double* ry = state->ry;
	double* rz = state->rz;
	double* vx = state->vx;
	double* vy = state->vy;
	double* vz = state->vz;
	double* energy = state->potentialEnergy;
	double* virial = state->virial;
	gid_type* label = state->label;
	SPECIES** species = state->species;
	GROUP** group = state->group;
	int nlocal = state->nlocal;
	gid_type nglobal = simulate->system->nglobal;

	PINFO_CODEC* pinfo = pinfoEncodeInit();
	unsigned pinfoFieldSize = bFieldSize(pinfoMaxIndex(pinfo));
	unsigned gidFieldSize = bFieldSize(mpiMaxVal(label, nlocal,COMM_LOCAL));

	mode =1; 
	switch (mode)
	{
	case 1:
		nfields = 11; 
		lrec = 9*4 + gidFieldSize + pinfoFieldSize;
		sprintf(field_names,"checksum id pinfo rx ry rz vx vy vz energy virial ");
		sprintf(field_types,"u4 b%1u b%1u f4 f4 f4 f4 f4 f4 f4 f4 ", gidFieldSize, pinfoFieldSize);
		break; 
	case 2: 
		nfields = 8; 
		lrec = 6*4 + gidFieldSize + pinfoFieldSize;
		sprintf(field_names,"checksum id pinfo rx ry rz energy virial ");
		sprintf(field_types,"u4 b%1u b%1u f4 f4 f4 f4 f4 ", gidFieldSize, pinfoFieldSize);
		break; 
	}
	if (Q!=NULL) 
	{
		sprintf(field_names+strlen(field_names),"Q%d C%d ",Lv[0],Lv[0]);
		sprintf(field_types+strlen(field_types),"f4 u4 ");
		lrec += 8;
		for (i=0;i<nL;i++) 
		{
			sprintf(field_names+strlen(field_names),"qn%d ",Lv[i]);
			sprintf(field_types+strlen(field_types),"f4 ");
			lrec += 4;
		}
		nfields +=(2+nL); 
	}
/* 	lrec =1;  */
/* 	while (lrec<4*nfields) lrec *= 2;  */
	unsigned char* line = (unsigned char*) ddcMalloc(lrec*sizeof(unsigned char));

	unsigned* offset = (unsigned*) ddcMalloc((nfields+1)*sizeof(unsigned));
	int nParse = parseFieldSizes(offset, field_types);
	assert(nParse == nfields);

	PioSet(file, "recordLength", lrec);
	PioSet(file, "datatype", FIXRECORDBINARY);
	PioSet(file, "numberRecords", nglobal);
	PioSet(file, "checksum", CRC32);
	PioSet(file, "nfields", nfields);
	PioSet(file, "field_names", field_names);
	PioSet(file, "field_types", field_types);
	PioSet(file, "misc_info", "groups =");
	for (unsigned ii=0; ii<pinfo->_nGroups; ++ii)
		PioSet(file, "misc_info", pinfo->_gNames[ii]);
	PioSet(file, "misc_info", ";\nspecies =");
	for (unsigned ii=0; ii<pinfo->_nSpecies; ++ii)
		PioSet(file, "misc_info", pinfo->_sNames[ii]);
	PioSet(file, "misc_info", ";\ntypes =");
	for (unsigned ii=0; ii<pinfo->_nTypes; ++ii)
		PioSet(file, "misc_info", pinfo->_tNames[ii]);
	PioSet(file, "misc_info", ";");
	      
	PioReserve(file, lrec*nlocal + 4096);
	if (getRank(0) == 0) write_fileheader(file, simulate, "bxyz" );
	double cLen = units_convert(1.0,NULL,"l"); 
	double time_convert = units_convert(1.0,NULL,"t"); 
	double energy_convert = units_convert(1.0,NULL,"eV"); 
	double cVel = cLen/time_convert; 
	for (i = 0; i < nlocal; i++)
	{
	   LONG64 pinfoIndex = pinfoEncode(group[i], species[i], pinfo);
	   float f4;
	   unsigned i4;
	   n=1;
	   bFieldPack(line+offset[n++], gidFieldSize,   label[i]);
	   bFieldPack(line+offset[n++], pinfoFieldSize, pinfoIndex);
	   double x = rx[i] ; double y = ry[i]; double z = rz[i];
	   backInBox(&x, &y, &z);
	   f4 = x*cLen;   copyBytes(line+offset[n++], &f4, 4);
	   f4 = y*cLen;   copyBytes(line+offset[n++], &f4, 4);
	   f4 = z*cLen;   copyBytes(line+offset[n++], &f4, 4);
	   if (mode == 1)
	   {
	      f4 = vx[i]*cVel; copyBytes(line+offset[n++], &f4, 4);
	      f4 = vy[i]*cVel; copyBytes(line+offset[n++], &f4, 4);
	      f4 = vz[i]*cVel; copyBytes(line+offset[n++], &f4, 4);
	   }
	   f4 = energy[i]*energy_convert;  copyBytes(line+offset[n++], &f4, 4);
	   f4 = virial[i]*energy_convert;  copyBytes(line+offset[n++], &f4, 4);
	   if (Q != NULL)
	   {
	      f4 = Q[i];  copyBytes(line+offset[n++], &f4, 4);
	      i4 = C[i];  copyBytes(line+offset[n++], &i4, 4);
	      for (k = 0; k < nL; k++)
	      {
		 f4 = qnorm[k][i];
		 copyBytes(line+offset[n++], &f4, 4);
	      }
	   }
	   i4 =  checksum_crc32_table(line+offset[1],lrec-offset[1]);
	   copyBytes(line, &i4, 4);
	   Pwrite(line, lrec, 1, (PFILE *)file);	   
	}
	ddcFree(offset);
	ddcFree(line);
	pinfoFree(pinfo);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
