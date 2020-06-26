#include "readCmds.h"

#include <string.h>
#include <unistd.h>
#include "ddcMalloc.h"
#include "mpiUtils.h"

#include "object.h"

static int isObject(char* line)
{
   for (unsigned ii=0; ii<strlen(line); ++ii)
      if (line[ii] == '{')
	 return 1;
   return 0;
}



int  readCMDS(char *filename) 
{
   FILE *file ;
   char line[256]; 
   int flag =0;
   line[0] = '\0';
   if (getRank(0) == 0) 
   {
      file = fopen(filename,"r");
      if (file)
      {
         fgets(line,255,file);
         if (! feof(file) && isObject(line))
         {
            object_compilefile(filename);
            flag |= NEW_OBJECT;
         }

         while (!feof(file))
         {
            size_t l=strlen(line); line[l-1] = '\0'; 
            printf("Received command \"%s\" from ddcMD_CMDS\n",line);
            if (strcmp(line,"checkpoint")==0) flag |= CHECKPOINT ; 
            if (strcmp(line,"kill")==0)       flag |= STOP ; 
            if (strcmp(line,"exit")==0)       flag |= STOP | CHECKPOINT;
            if (strcmp(line,"profile")==0)    flag |= DUMP_PROFILE;
            if (strcmp(line,"hpm")==0)        flag |= HPM_PRINT;
            if (strcmp(line,"analysis")==0)   flag |= DO_ANALYSIS;

            fgets(line,255,file);
         }
         fclose(file);
         truncate(filename, 0);
      }
   }
   MPI_Bcast(&flag,1,MPI_INT,0, COMM_LOCAL); 

   return flag; 
}



/** Rescans a variety of parameters out of the object database.  This
 * function is intended to be called when an obect has been read from
 * ddcMD_CMDS. */
void object_rescan(DDC* ddc, SIMULATE* simulate)
{
   int nelements;
   OBJECT* obj;
   double d[1024];
   char*  ignore;
   ignore = (char *) IGNORE_IF_NOT_FOUND;

   obj = object_find("ddc", "DDC");
   if (obj)
   {
      nelements = object_get(obj, "dx", d, DOUBLE, 1024, "-1.0"); 
      if (d[0] > 0.0) {ddc->lx = nelements; ddc->dx = ddcRealloc(ddc->dx, ddc->lx*sizeof(double)); for (int i=0; i<ddc->lx; i++) ddc->dx[i] = d[i];}
      nelements = object_get(obj, "dy", d, DOUBLE, 1024, "-1.0"); 
      if (d[0] > 0.0) {ddc->ly = nelements; ddc->dy = ddcRealloc(ddc->dy, ddc->ly*sizeof(double)); for (int i=0; i<ddc->ly; i++) ddc->dy[i] = d[i];}
      nelements = object_get(obj, "dz", d, DOUBLE, 1024, "-1.0"); 
      if (d[0] > 0.0) {ddc->lz = nelements; ddc->dz = ddcRealloc(ddc->dz, ddc->lz*sizeof(double)); for (int i=0; i<ddc->lz; i++) ddc->dz[i] = d[i];}
      nelements = object_get(obj, "updateRate", &ddc->updateRate, INT, 1, ignore);
   }

   obj = object_find("test", "SIMULATE");
   if (obj)
   {
      object_get(obj, "checkpointrate", &simulate->checkpointrate, INT, 1, ignore);
      object_get(obj, "snapshotrate",   &simulate->snapshotrate,   INT, 1, ignore);
      object_get(obj, "maxloop",        &simulate->maxloop,        U64, 1, ignore);
   }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
