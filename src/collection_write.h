#ifndef COLLECTION_WRITE_HH
#define COLLECTION_WRITE_HH

#include "collection.h"
#include "simulate.h"
#include "pio.h"

void collection_writeMode(int);
void collection_nCopies(int n);


void collection_writeBLOCK(SIMULATE*simulate, PFILE* file);
void collection_writeBLOCK_binary(SIMULATE*simulate, PFILE* file);
void collection_writeBXYZ(SIMULATE*simulate, PFILE* file);

#endif // #ifndef COLLECTION_WRITE_HH


/* Local Variables: */
/* tab-width: 3 */
/* End: */
