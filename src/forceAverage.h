#include "object.h"

typedef struct forceAverage_parms_st
{
   
   // standard analysis parameters
   int       _eval_rate, _outputrate;
   int       _nFiles;
   char*     _filename;
   char*     _data_type;
   char*     _parms_info;

   // Record forces for these particles and their gid.
   int        _nParticles;
   gid_type * _idParticle;

   FILE*     _file;
   THREE_VECTOR * _forbar;
   int _nSnapSum;

} FORCE_AVERAGE_PARMS;
