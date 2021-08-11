#include <string.h>
#include <stdlib.h>

#include "object.h"
#include "ddcMalloc.h"
#include "error.h"
#include "units.h"
#include "utilities.h"
#include "ioUtils.h"
#include "object.h"
#include "mpiUtils.h"
#include "masters.h"
#include "commandLineOptions.h"
void objectSetup(void *parms, MPI_Comm comm)
{
   COMMONPARMS commonParms =((DEFAULTMASTERPARMS *)parms)->common;
   
   timestamp("Start objectSetup");
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   if (myRank == 0)
   {
      char *objectFiles = strdup(commonParms.objectFile); 
      char *objectFile = strtok(objectFiles," "); 
      while(objectFile != NULL) 
      {
         if (filetest(objectFile, S_IFREG) != 0)
         {
            printf("objectfile=%s does not exist or wrong type\n", objectFile);
            abortAll(1);
         }
         objectFile = strtok(NULL," "); 
      }
      free(objectFiles); 
      if (filetest(commonParms.restartFile, S_IFREG) != 0)
      {
         printf("restart=%s does not exist or wrong type\n", commonParms.restartFile);
         abortAll(1);
      }
      char* string = (char *)ddcMalloc(strlen(commonParms.objectFile) + 1 + strlen(commonParms.restartFile) + 1+ strlen("potential.data")+2);
      sprintf(string, "%s %s", commonParms.objectFile, commonParms.restartFile);
      struct stat buf;
      if (stat("potential.data", &buf ) == 0) strcat(string, " potential.data");
      object_set("files", string);
      ddcFree(string);

      char filename[1024];
      OBJECT *c,*obj;
      object_compilestring("routineManager ROUTINEMANAGER { routines=simulate SIMULATE; layouts=-1;} " ) ;
      object_compilestring("ddc DDC {} " ) ;
      object_compilestring("energyInfo ENERGYINFO {} " ) ;
      object_compilestring("loadbalance LOADBALANCE {}");
      object_compile();
      c = object_find("collection", "COLLECTION");
      int nelements = object_get(c, "files", &string, STRING, 1, " ");
      if (nelements > 0) 
      {
         strcpy(filename,string);
         string = strchr(filename,'#');
         if (string) 
         {
            strcpy(string+1,"000000");
            object_compilefilesubset(filename,0,0);
            obj = object_find2("particle", "FILEHEADER",IGNORE_IF_NOT_FOUND);
            if (!obj) obj = object_find2("bxyz", "FILEHEADER",IGNORE_IF_NOT_FOUND);
            if (!obj) obj = object_find2("header", "FILEHEADER",IGNORE_IF_NOT_FOUND);
            if (obj) 
            {
               obj->name="file";
               object_get(c, "headerLength", &string, STRING, 1, "0");
               object_replacekeyword(obj,"headerLength",string);
            }
         }
      }
   }
   object_Bcast(0, comm);
   WAIT(-1); 
   timestamp("Finish objectSetup");
}
