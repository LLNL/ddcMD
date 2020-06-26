#include "commandLineOptions.h"
#include "ddcMalloc.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "utilities.h"

/*
void commonDefault(void **opt_, char *masterName)
{
   COMMONOPT *opt = (COMMONOPT *)(*opt_); 
   opt->masterName = strdup(masterName); 
	opt->objectFile = "object.data";
	opt->restartFile = "restart";
	opt->simulateName = NULL;
}
void commonOpt(char *option, char *value, void** opt_)
{
      COMMONOPT *opt = (COMMONOPT *)(*opt_); 
		if (strcmp(option, "-s") == 0) opt->simulateName = strdup(value);
		if (strcmp(option, "-r") == 0) opt->restartFile  = strdup(value);
		if (strcmp(option, "-o") == 0) opt->objectFile   = strdup(value);
}
void simulateOpt(char *option, char *value, void ** opt_)
{
   char *masterName = "simulateMaster";
   SIMULATEOPT *opt = (SIMULATEOPT *)(*opt_); 
   if (opt == NULL)  
   {
      opt = (SIMULATEOPT *)ddcMalloc(sizeof(SIMULATEOPT)); 
      commonDefault(opt_,masterName);
   }
   assert(strcmp(opt->masterName,masterName) == 0);
   commonOpt(option, value, opt_);
   if (strcmp(option, "-STOP_TIME") == 0) opt->stopTime=strdup(value); 
}

void eightFoldOpt(char *option, char *value, void** opt_)
{
   char *masterName = "eightFoldMaster";
   EIGHTFOLDOPT *opt = (EIGHTFOLDOPT *)(*opt_); 
   if (opt == NULL)  
   {
      opt = (EIGHTFOLDOPT *)ddcMalloc(sizeof(EIGHTFOLDOPT)); 
      commonDefault(opt_,masterName);
   }
   assert(strcmp(opt->masterName,masterName) == 0);
   commonOpt(option, value, opt_);
}
void readWriteOpt(char *option, char *value, void** opt_)
{
   char *masterName = "readWriteMaster";
   READWRITEOPT *opt = (READWRITEOPT *)(*opt_); 
   if (opt == NULL)  
   {
      opt = (READWRITEOPT *)ddcMalloc(sizeof(READWRITEOPT)); 
      commonDefault(opt_,masterName);
   }
   assert(strcmp(opt->masterName,masterName) == 0);
   commonOpt(option, value, opt_);
}
*/
/*
	opt.stopTime = NULL;
*/

COMMAND_LINE_OPTIONS  parseCommandLine(int argc, char** argv)
{
   COMMAND_LINE_OPTIONS  opt; 
   int nOptions =0; 
	for (int i=1; i<argc;i++) if (argv[i][0] == '-') nOptions++; 
   opt.options=ddcMalloc(nOptions*sizeof(char *));
   opt.values=ddcMalloc(nOptions*sizeof(char *));
   opt.masterName = strdup("simulateMaster"); 
   int i=1; 
   if (argc > 1) 
   {
   if (argv[i][0] != '-' )  
   {
      free(opt.masterName); 
      opt.masterName = strdup(argv[1]); 
      astrcat(opt.masterName,"Master"); 
      i++;
   }
   }
   opt.nOptions = nOptions; 
   int n =0; 
	for (; i<argc;)
	{
      assert(argv[i][0] == '-');
      
      char *option = strdup(argv[i]); 
      char *value  = strdup(""); 
      int j=i+1; 
      while (j<argc && argv[j][0] != '-' ) 
      {
         if (value[0] != '\0') value = astrcat(value," "); 
         value = astrcat(value,argv[j++]); 
      }
      opt.options[n] = strdup(option); 
      opt.values[n] = strdup(value); 
      n++; 
      


/*
	if (!strcmp(option, "transform")       {opt.masterName  = strdup("transformMaster");       opt.masterValue  = strdup(value); }
		if (!strcmp(option, "analysis")        {opt.masterName  = strdup("analysisMaster");        opt.masterValue  = strdup(value); }
		if (!strcmp(option, "thermalize")      {opt.masterName  = strdup("thermalizeMaster");      opt.masterValue  = strdup(value); }
		if (!strcmp(option, "integrationTest") {opt.masterName  = strdup("integrationTestMaster"); opt.masterValue  = strdup(value); }
      if (!strcmp(option, "unitTest")        {opt.masterName  = strdup("unitTestMaster");        opt.masterValue  = strdup(value); }
*/
      free(option); 
      free(value); 
      i=j; 
	}
	return opt;
}
/* Local Variables: */
/* tab-width: 3 */
/* End: */
