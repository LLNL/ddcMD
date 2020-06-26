#include "io.h"
#include "sCatPrintf.h"

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define  _LARGE_FILE
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <libgen.h> // dirname
#include "three_algebra.h"
#include "ddcMalloc.h"
#include "units.h"
#include "collection_write.h"
#include "ioUtils.h"
#include "mpiUtils.h"
#include "utilities.h"
#include "integrator.h"
#include "object.h"
#include "random.h"
#include "format.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* PATH_MAX is not defined in limits.h on some platforms */

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

enum CHECKPOINTUNITS_ENUM
{
   LENGTH, MASS, TIME, CURRENT, TEMPERATURE, AMOUNT, LUMINOUS_INTENSITY
};
static char *checkpointUnitNames[7];
char *version_get(char *);

void ddc_writePXYZ(void* ddc, SIMULATE* simulate, PFILE* file);
void checkpointUnits(char *lengthUnit, char *massUnit, char *timeUnit, char *currentUnit, char* temperatureUnit, char *amountUnit, char *luminous_intensityUnit)
{
   checkpointUnitNames[LENGTH]=strdup(lengthUnit); 
   checkpointUnitNames[MASS]=strdup(massUnit); 
   checkpointUnitNames[TIME]=strdup(timeUnit); 
   checkpointUnitNames[CURRENT]=strdup(currentUnit); 
   checkpointUnitNames[TEMPERATURE]=strdup(temperatureUnit); 
   checkpointUnitNames[AMOUNT]=strdup(amountUnit); 
   checkpointUnitNames[LUMINOUS_INTENSITY]=strdup(luminous_intensityUnit); 
}
void writeRestart(SIMULATE*simulate, int restartLink)
{
   PFILE *pfile;
   FILE *file;
	char filename[512];
	int rc, id, i, sizegroup, ngroup;
	timestamp("Start Writing Restart");
	id = getRank(0);
	sprintf(filename, "%s/atoms", simulate->snapshotdir);
	pfile = Popen(filename, "w", COMM_LOCAL);
	if (simulate->checkpointmode == BINARY_CHECKPOINT)
	   collection_writeBLOCK_binary(simulate, pfile);
	else
	   collection_writeBLOCK(simulate, pfile);
	Pget(pfile, "sizegroup", &sizegroup);
	Pget(pfile, "ngroup", &ngroup);
	rc = Pclose(pfile);
	if (id == 0)
	{
		sprintf(filename, "%s/restart", simulate->snapshotdir);
		file = fopen(filename, "w");
		fprintf(file, "%s SIMULATE { run_id=0x%08x; loop=%"PRId64"; time=%f fs;}\n", simulate->name, simulate->run_id, simulate->loop, 
         units_convert(simulate->time,NULL,"t"));fflush(file); 
		box_write(simulate->system->box, file);
		simulate->integrator->writedynamic(simulate->integrator, file); fflush(file); 
		for (i = 0; i < simulate->system->ngroup; i++) 
		{
			GROUP *g ;
			g = simulate->system->group[i];
			if (g->write_dynamics != NULL) g->write_dynamics(g,file); 
			fflush(file); 
		}
		for (i = 0; i < simulate->system->npotential; i++) 
		{
			POTENTIAL *p ;
			p = simulate->system->potential[i];
			if (p->write_dynamics != NULL) p->write_dynamics(p,file); fflush(file); 
		}
      int nDynamic = simulate->dynamicWrite.n; 
      DYNAMICRESTARTINFO *info  = simulate->dynamicWrite.info; 
      for (i =0; i < nDynamic; i++) 
      {
         info[i].fcn(info[i].parms,file);
      }

		fprintf(file, "%s COLLECTION { size=%"PRIu64"; files=%s/atoms#;}\n",
			simulate->system->collection->name, simulate->system->nglobal, simulate->snapshotdir); fflush(file); 
		fclose(file);
		if (rc == 0 && restartLink) 
		{
			unlink("restart");
			symlink(filename,"restart");
		}
	}
	if (getSize(0) != 1) writePXYZ(simulate->ddc, simulate);
	timestamp("Finished Writing Restart");
}
char *CreateSnapshotdir(SIMULATE *simulate, char* dirname)
{
	int id; 
	id = getRank(0);
	if (dirname == NULL)
	{
	   int nameLength =
		  strlen(simulate->atomsdir) /* "%s" for simulate->atomsdir */
		  + 9 /* "/snapshot" */
		  + 1 /* "." */
		  + MAX(16,(log(fabs(simulate->maxloop*1.0))/log(10.0)+1  )) /* "loopFormat" */
		  + 1 /* Terminating zero byte ('\0') */;
		simulate->snapshotdir = ddcRealloc(simulate->snapshotdir,nameLength); 
	   sprintf(simulate->snapshotdir, "%s/snapshot.", simulate->atomsdir);
	   sprintf(simulate->snapshotdir+strlen(simulate->snapshotdir), loopFormat(), simulate->loop);
	}
	else
	{
	   int nameLength = strlen(simulate->atomsdir) +  1 + strlen(dirname)   + 1; 
		simulate->snapshotdir = ddcRealloc(simulate->snapshotdir,nameLength); 
	   sprintf(simulate->snapshotdir, "%s/%s", simulate->atomsdir,dirname);
	}
	
	if (id == 0 ) DirTestCreate(simulate->snapshotdir);
   WAIT(0);
   usleep(10000);
   return simulate->snapshotdir; 
}

void writeBXYZ(SIMULATE*simulate)
{
   PFILE *file;
   char filename[512];
   timestamp("Start Writing bxyz");
   sprintf(filename, "%s/bxyz", simulate->snapshotdir);
   file = Popen(filename, "w", COMM_LOCAL);
   collection_writeBXYZ(simulate, file);
   Pclose(file);
   WAIT(0);
   timestamp("Finished Writing bxyz");
}

void writePXYZ(DDC* ddc, SIMULATE* simulate)
{
   PFILE *file;
   char filename[512];
   timestamp("Start Writing pxyz");
   sprintf(filename, "%s/pxyz", simulate->snapshotdir);
   file = Popen(filename, "w", COMM_LOCAL);
   ddc_writePXYZ(ddc, simulate, file);
   Pclose(file);
   WAIT(0);
   timestamp("Finished Writing pxyz");
}

void writeRestart8(SIMULATE*simulate)
{
   PFILE *file;
   char filename[512], dirname[512];
   int id, i, j,sizegroup, ngroup, recordLength;
   STATE *state;
   gid_type *label_new;
   double *rx, *ry, *rz, *vx, *vy, *vz, volume;
   gid_type shift[] = { 1,5,2,7,3,0,4,6}; 
   gid_type *label;
   char *field = NULL;
   COLLECTION *c;
   void (*blockWriteFunction) (SIMULATE*simulate, PFILE *file);
   char mode[32];
   if (simulate->checkpointmode == BINARY_CHECKPOINT)
   {
      blockWriteFunction = collection_writeBLOCK_binary;
      sprintf(mode, "%s", "FIXRECORDBINARY");
   }
   else
   {
      blockWriteFunction = collection_writeBLOCK;
      sprintf(mode, "%s", "FIXRECORDASCII");
   }

   timestamp("Start Writing Restart");
   id = getRank(0);
   strcpy(dirname, "snapshot.initial");
   if (id == 0) DirTestCreate(dirname);
   WAIT(0);
   sprintf(filename, "%s/atoms", dirname);
   file = Popen(filename, "w", COMM_LOCAL);
   THREE_MATRIX h0 = box_get_h(NULL);
   double xSide = h0.xx;
   double ySide = h0.yy;
   double zSide = h0.zz;
   if ( h0.xy!=0 || h0.xz!=0 || h0.yx!=0 || h0.yz!=0 || h0.zx!=0 || h0.zy!=0)
   {
      if ( getRank(0) == 0)
         printf("8fold supports only orthorhombic boxes.\n");
      exit(3);
   }
   state = simulate->system->collection->state;
   c = simulate->system->collection;
   rx = state->rx;
   ry = state->ry;
   rz = state->rz;
   vx = state->vx;
   vy = state->vy;
   vz = state->vz;
   label = state->label;
   box_get(simulate->system->box, VOLUME, (void *)&volume);
   volume *= 8.0;
   box_put(simulate->system->box, VOLUME, (void *)&volume);
   collection_writeMode(WHEADER8FOLD);
   label_new = ddcMalloc(c->size*sizeof(*label_new));
   for (i = 0; i < c->size; i++) label_new[i]  = 8*label[i]; 
   LONG64 seed = generateRandomSeed();
   double delta=0.0; 
   for (i = 0; i < c->size; i++)
   {
      j = i % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rx[i] -= xSide*0.5;
      ry[i] -= ySide*0.5;
      rz[i] -= zSide*0.5;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 000 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+1) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rx[i] += xSide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   collection_writeMode(NOWHEADER);
   blockWriteFunction(simulate, file);	/* 100 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+2) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      ry[i] += ySide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 110 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+3) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rz[i] += zSide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 111 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+4) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rx[i] -= xSide;
      rz[i] -= zSide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 010 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+5) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rz[i] += zSide;
      ry[i] -= ySide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 001 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+6) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rx[i] += xSide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 101 */
   for (i = 0; i < c->size; i++)
   {
      j = (i+7) % 8  ; 
      label[i] =  label_new[i] +shift[j]; 
      PRAND48_STATE handle = prand48_init(label[i], seed, 0x96ab45cd12ef76llu);
      rx[i] -= xSide;
      ry[i] += ySide;
      vx[i] += delta*gasdev0(&handle);
      vy[i] += delta*gasdev0(&handle);
      vz[i] += delta*gasdev0(&handle);
      if (state->group[i]->parse != NULL) state->group[i]->parse(state->group[i], field, i);
   }
   blockWriteFunction(simulate, file);	/* 011 */
   Pget(file, "sizegroup", &sizegroup);
   Pget(file, "ngroup", &ngroup);
   Pget(file, "recordLength", &recordLength);
   Pclose(file);
   if (id == 0)
   {
      FILE* rfile;
      sprintf(filename, "%s/restart", dirname);
      rfile = fopen(filename, "w");
      fprintf(rfile, "%s SIMULATE { run_id=0x%08x; loop=%"PRId64"; time=%.4e fs;}\n", simulate->name, simulate->run_id, simulate->loop, units_convert(simulate->time,NULL,"t"));

      box_write(simulate->system->box, rfile);
      fprintf(rfile, "%s COLLECTION { mode=%s; recordLength=%d; size=%"PRIu64"; files=%s/atoms#;}\n",
            simulate->system->collection->name, mode, recordLength, simulate->system->nglobal*8, dirname);
      fclose(rfile);
   }
   timestamp("Finish Writing Restart");
   collection_writeMode(WHEADER);
   ddcFree(label_new);
}

void write_fileheader(PFILE *file, SIMULATE *simulate, char *name )
{
//#define HEADERLENGTH  8192 
   double *hptr;
   int nfiles,nfields,lrec,key;
   gid_type nrecord;
   //char line[HEADERLENGTH],*checksumName,*datatypeName,*field_names,*field_types, *field_units, *field_format, *misc_info,*create_time; 
   char *checksumName,*datatypeName,*field_names,*field_types, *field_units, *field_format, *misc_info,*create_time; 
   char *line=NULL;
   memcpy(&key, "1234", 4);
   box_get(NULL, HO_PTR, (void *)&hptr);
   THREE_VECTOR reducedcorner = box_get_reducedcorner(NULL); 
   Pget(file,"ngroup",&nfiles);
   Pget(file,"numberRecords",&nrecord);
   Pget(file,"recordLength",&lrec);
   Pget(file,"datatypeName",&datatypeName);
   Pget(file,"checksumName",&checksumName);
   Pget(file,"nfields",&nfields);
   Pget(file,"field_names",&field_names);
   Pget(file,"field_types",&field_types);
   Pget(file,"field_units",&field_units);
   Pget(file,"field_format",&field_format);
   Pget(file,"misc_info",&misc_info);
   create_time=timestamp_string(); 

   sCatPrintf(&line, "%s FILEHEADER {type=MULTILINE; datatype=%s; checksum=%s; create_time=%s; run_id=0x%08x;\n",
         name,datatypeName,checksumName,create_time,simulate->run_id);
   sCatPrintf(&line, "code_version=%s; srcpath=%s;\n", version_get("version"),version_get("srcpath"));
   sCatPrintf(&line, "loop=%"PRId64"; time=%f fs;\n", simulate->loop,units_convert(simulate->time,NULL,"t"));
   sCatPrintf(&line, "nfiles=%d; nrecord=%"PRIu64"; lrec=%d; nfields=%d; endian_key=%d;\n", nfiles,nrecord, lrec, nfields, key);
   sCatPrintf(&line, "field_names=%s;\n", field_names);
   sCatPrintf(&line, "field_types=%s;\n", field_types);
   if (field_units != NULL ) sCatPrintf(&line, "field_units=%s;\n", field_units);
   if (field_format != NULL ) sCatPrintf(&line, "field_format=%s;\n", field_format);
   sCatPrintf(&line, "reducedcorner=%21.14f %21.14f %21.14f;\n", reducedcorner.x,reducedcorner.y,reducedcorner.z);
   double length_convert = units_convert(1.0,NULL,"l"); 
   sCatPrintf(&line, "h=%21.14f %21.14f %21.14f\n", hptr[0]*length_convert,hptr[1]*length_convert,hptr[2]*length_convert);
   sCatPrintf(&line, "  %21.14f %21.14f %21.14f\n", hptr[3]*length_convert,hptr[4]*length_convert,hptr[5]*length_convert);
   sCatPrintf(&line, "  %21.14f %21.14f %21.14f Ang;\n", hptr[6]*length_convert,hptr[7]*length_convert,hptr[8]*length_convert);
   if(strlen(misc_info) > 0) sCatPrintf(&line, "%s\n",misc_info);
   sCatPrintf(&line, "}\n");
   int headerlength=sCatPrintf(&line, " \n\n");
   /*
      for (i=0;i<simulate->system->ngroup;i++)
      {

      GROUP *gp = simulate->system->group[i]; 
      sCatPrintf(&line, "%s GROUP {%s}\n",gp->name,gp->value); 
      }
      */
   Pwrite(line,headerlength,1,(PFILE *)file);
   free(line);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
