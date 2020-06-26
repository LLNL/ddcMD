#include "readPXYZ.h"

#include <mpi.h>
#include <string.h>
#include <assert.h>
#include <libgen.h> // dirname

#include "utilities.h"
#include "mpiUtils.h"
#include "pio.h"
#include "units.h"
#include "heap.h"
#include "ddcMalloc.h"
#include "ioUtils.h"
#include "mpiTypes.h"


static int ddc_readPXYZ(DDC* ddc, PFILE* file);


/**
 *  Reads a pxyz file and stores the domain center data in ddc->domain.
 *
 *  Returns zero on successful read, non-zero otherwise.  You should
 *  check the return value from this function and provide alternate
 *  initialization of the domain centers if this read fails.  Possible
 *  failure modes include restarting the job on a different number of
 *  tasks from the previous run.
 *
 *  In this routine we make as few assumptions about which objects have
 *  been initialized as possible.  In particular, we do not assume that
 *  either the SYSTEM or the COLLECTION have been initialized.  Instead
 *  we chain our way down the object database hierachy to get the name
 *  of the directory from which the collection data will be read (if it
 *  is to be read).
 *
 *  We *do* assume that the DDC has been sufficiently initialized that
 *  parent is set and memory has been allocated for the domain array.
 *  
 */
int readPXYZ(DDC* ddc)
{
   timestamp("Start reading pxyz");
	OBJECT* parentObj = ddc->parent;
	char* sysName;
	object_get(parentObj, "system", &sysName, STRING, 1, "");
	OBJECT* sysObj = object_find(sysName, "SYSTEM");
	char* colName;
	object_get(sysObj, "collection", &colName, STRING, 1, "");
	OBJECT* colObj = object_find(colName, "COLLECTION");
	char* atomsName;
	if (object_testforkeyword(colObj, "files") == 0)
		 return -1;
	object_get(colObj, "files", &atomsName, STRING, 1, "");
	
	char* tmpName = strdup(atomsName);
	char* dirName = dirname(tmpName);
	int len = strlen(dirName);
	char filename[len+7];
	sprintf(filename, "%s/pxyz#", dirName);
	char masterFile[len+13];
	sprintf(masterFile, "%s/pxyz#000000", dirName);
	int pxyzExists;
	if (getRank(0) == 0)
		pxyzExists = filetest(masterFile, S_IFREG);
	MPI_Bcast(&pxyzExists, 1, MPI_INT, 0, COMM_LOCAL);
	if (pxyzExists != 0)
		return -2;

	
	PFILE* file = Popen(filename, "r", COMM_LOCAL);
	assert(file != 0);
	ddcFree(tmpName);
	
	// let someone else do the heavy lifting.
	int status = ddc_readPXYZ(ddc, file);
   timestamp("Finished reading pxyz");
	return status;
}


/**
 *  Reads the positions of the domain centers from the supplied PFILE
 *  and populates the centers in ddc->domain (radius information is not
 *  saved in pxyz so it is simply set to zero by this routine).  This
 *  routine does broadcast the locations of all centers to all tasks.
 *
 *  This routine will only populate ddc->domain if it is running on the
 *  same number of tasks as when the pxyz file was written.  If a
 *  different number of tasks is detected this routine will return the
 *  number of centers that are written in the pxyz file.
 *
 *  You might think that since we only run this code when we are reading
 *  on the same number of tasks as the pxyz files were written on that
 *  we are guaranteed that each pio will deliver exactly one record to
 *  each task to process.  However, this is not the case.  When the
 *  number of tasks is not exactly divisible by the number of files the
 *  logic for splitting up the excess tasks is different in the reading
 *  and writing cases.  Hence, this code must handle the case that some
 *  tasks get multiple records while some get none.
 *
 *  Returns zero on a successful read, non-zero otherwise.
 */
static int ddc_readPXYZ(DDC* ddc, PFILE* file)
{
   int nTasks;
   MPI_Comm_size(file->comm, &nTasks);
   int myRank;
   MPI_Comm_rank(file->comm, &myRank);
   assert(myRank == file->id);
   assert(ddc->domains.size == nTasks);
   assert(file->datatype == FIXRECORDASCII);
   
   if ( (int)file->numberRecords != ddc->domains.size)
      return file->numberRecords;

   unsigned nRecords = file->bufsize / file->recordLength;
   assert(file->bufsize % file->recordLength == 0);


   char** fieldUnits = NULL;
   unsigned nUnits = object_getv(
      file->headerObject, "field_units", (void*)&fieldUnits,
      STRING, ABORT_IF_NOT_FOUND);
   assert(nUnits >= 4);
   
   DOMAINX dIn[nRecords];
   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      char buf[file->recordLength];
      buf[0] = '\0';
      Pfgets(buf, file->recordLength, file);
      if (strlen(buf) < 60)
      { // deal with old, broken pxyz files with extra newlines.
	 buf[0] = '\0';
	 Pfgets(buf, file->recordLength, file);
      }
      
      unsigned id;
      THREE_VECTOR r;
      sscanf(buf, "%u %lf %lf %lf", &id, &r.x, &r.y, &r.z);
      //printf("ddt %u: %u %lf %lf %lf\n", myRank, id, r.x, r.y, r.z);

      r.x *= units_convert(1, fieldUnits[1], NULL);
      r.y *= units_convert(1, fieldUnits[2], NULL);
      r.z *= units_convert(1, fieldUnits[3], NULL);
      //printf("%u %lf %lf %lf\n", id, r.x, r.y, r.z);
      dIn[ii].radius = 0;
      dIn[ii].center = r;
   }

   unsigned recvCntBlk = 0;
   unsigned displsBlk = 0;
   int* recvCnt = NULL;
   int* displs = NULL;
   if (myRank == 0)
   {
      recvCnt = heapGet(&recvCntBlk);
      heapEndBlock(recvCntBlk, nTasks*sizeof(unsigned));
      displs = heapGet(&displsBlk);
      heapEndBlock(displsBlk, nTasks*sizeof(unsigned));
   }
   MPI_Gather(&nRecords, 1, MPI_INT, recvCnt, 1, MPI_INT, 0, file->comm);
   if (myRank == 0)
   {
      displs[0] = 0;
      for (int ii=0; ii<nTasks-1; ++ii)
	 displs[ii+1] = displs[ii] + recvCnt[ii];
   }
   
   MPI_Gatherv(dIn, nRecords, domainx_MPIType(),
	       ddc->domains.domains, recvCnt, displs, domainx_MPIType(),
	       0, file->comm);
   MPI_Bcast(ddc->domains.domains, nTasks, domainx_MPIType(), 0, file->comm);

   if (myRank == 0)
   {
      heapFree(displsBlk);
      heapFree(recvCntBlk);
   }
   for (unsigned ii=0; ii<nUnits; ++ii)
      ddcFree(fieldUnits[ii]);
   ddcFree(fieldUnits);

/*    for (int ii=0; ii<nTasks; ++ii) */
/*       fprintf(ddcfile, "%u: (%f, %f, %f)  %f\n", ii, */
/* 	      ddc->domains.domains[ii].center.x, */
/* 	      ddc->domains.domains[ii].center.y, */
/* 	      ddc->domains.domains[ii].center.z, */
/* 	      ddc->domains.domains[ii].radius); */
/*    fflush(ddcfile); */

   
   return 0;
}
