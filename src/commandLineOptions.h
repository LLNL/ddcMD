#ifndef COMMANDLINEOPTIONS_H
#define COMMANDLINEOPTIONS_H

#include <time.h>
#include <mpi.h>

/** Process and store command line options. */


/*
typedef struct CommandLineOptions_st
{
   double thermalTemp;

   char*  masterName; 
   char*  masterValue; 
   char*  simulateName;
   char*  restartFile;
   char*  objectFile;
   char *stopTime;
   char *value; 
   
} COMMAND_LINE_OPTIONS ;
*/
typedef struct CommandLineOptions_st 
{
   char *masterName;
   int nOptions;
   char **options; 
   char **values; 
} COMMAND_LINE_OPTIONS ;

COMMAND_LINE_OPTIONS  parseCommandLine(int argc, char** argv);

#endif // #ifndef COMMAND_LINE_OPTIONS 
