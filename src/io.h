#ifndef IO_H
#define IO_H

#include "simulate.h"
#include "ddc.h"
#include "pio.h"



char *CreateSnapshotdir(SIMULATE *simulate, char* dirname);
void writeRestart(SIMULATE*simulate,int restartLink);
void writeRestart8(SIMULATE*simulate);
void writeBXYZ(SIMULATE*simulate);
void writePXYZ(DDC* ddc, SIMULATE* simulate);

void write_fileheader(PFILE* file, SIMULATE* simulate, char* name );
void checkpointUnits(char *lengthUnit, char *massUnit, char *timeUnit, char *currentUnit, char* temperatureUnit, char *amountUnit, char *luminous_intensityUnit);



#endif // #ifndef IO_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
