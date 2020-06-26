#include <string.h>
#include <stdlib.h>

#include "utilities.h"
#include "units.h"
#include "commandLineOptions.h"
#include "ddcMalloc.h"
#include "masters.h"
//

int getRank(int); 
COMMONPARMS commonFactory(COMMAND_LINE_OPTIONS opt)
{
   COMMONPARMS  parms={"simulate","object.data","restart"}; 
   for (int i=0;i<opt.nOptions;i++)
   {
      if (strcmp(opt.options[i],"-s")==0) parms.simulateName = strdup(opt.values[i]);
      if (strcmp(opt.options[i],"-o")==0) parms.objectFile =   strdup(opt.values[i]);
      if (strcmp(opt.options[i],"-r")==0) parms.restartFile =  strdup(opt.values[i]);
   }
   return parms; 
}
MASTER masterFactory(COMMAND_LINE_OPTIONS opt)
{
   MASTER master={NULL,NULL}; 
   COMMONPARMS commonParms = commonFactory(opt); 
   if (strcmp(opt.masterName,"simulateMaster")==0) 
   {
      SIMULATEMASTERPARMS *parms =  (SIMULATEMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      parms->stopTime = -1; 
      for (int i=0;i<opt.nOptions;i++)
         if (strcmp(opt.options[i],"-STOP_TIME")==0) parms->stopTime = convert_timestring(opt.values[i]);

      master.fcn=simulateMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"analysisMaster")==0) 
   {
      ANALYSISMASTERPARMS *parms = (ANALYSISMASTERPARMS *)ddcMalloc(sizeof(ANALYSISMASTERPARMS)); 
      parms->common = commonParms; 
      parms->energyCalc = 1; 
      for (int i=0;i<opt.nOptions;i++)
         if (strcmp(opt.options[i],"-noEnergy")==0) parms->energyCalc = 0;

      master.fcn=simulateMaster;
      master.fcn=analysisMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"testPressureMaster")==0) 
   {
      DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      master.fcn=testPressureMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"testForceMaster")==0) 
   {
      DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common=commonParms;
      master.fcn=testForceMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"eightFoldMaster")==0) 
   {
      DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      master.fcn=eightFoldMaster;
      parms->common = commonParms; 
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"readWriteMaster")==0) 
   {
      DEFAULTMASTERPARMS *parms = master.parms = (DEFAULTMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      master.fcn=readWriteMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"thermalizeMaster")==0) 
   {
      THERMALIZEMASTERPARMS *parms = (THERMALIZEMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
	   parms->thermalTemp = -1.0;
      for (int i=0;i<opt.nOptions;i++)
      {
         if (strcmp(opt.options[i],"-T")==0) 
         {
            char *units=NULL; 
            double thermalTemp = strtod(opt.values[i],&units);
            if (units[0]=='\0') units="T"; 
            parms->thermalTemp = units_convert(thermalTemp,units,NULL);
            break; 
         }
      }
      master.fcn=thermalizeMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"transformMaster")==0) 
   {
      DEFAULTMASTERPARMS *parms = (DEFAULTMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      master.fcn=transformMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"integrationTestMaster")==0) 
   {   
      SIMULATEMASTERPARMS *parms =  (SIMULATEMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      parms->stopTime = -1; 

      master.fcn=integrationTestMaster;
      master.parms=parms; 
   }
   if (strcmp(opt.masterName,"unitTestMaster")==0) 
   {   
      SIMULATEMASTERPARMS *parms =  (SIMULATEMASTERPARMS *)ddcMalloc(sizeof(SIMULATEMASTERPARMS)); 
      parms->common = commonParms; 
      parms->stopTime = -1; 

      master.fcn=unitTestMaster;
      master.parms=parms; 
   }
   return master; 
}
