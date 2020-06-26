#ifndef TRANSECTMORPH_H
#define TRANSECTMORPH_H
#include "group.h"
#include "box.h"
#include "transform.h"

enum NORMALTYPE {IN_UNITS,IN_LATTICE};

typedef struct transect_morph_st
{
 // The set of all planes is defined by a normal and
 // a list of scalar positions along that normal.
 // The normal is parallel to the index'th cell axis.
 // "Before" planes are to be shifted to "After" positions.
 // Planes are stored in positional order, both Before and After.
 // Before planes are assumed to all lie inside the unit cell.
 // Wraparound of PBC is handled by letting After positions
 // lie outside the unit cell.
 // After planes are not allowed to pass through one another;
 // PBC complicate this constraint slightly...

 int index;
 enum NORMALTYPE units;
 double *positionBefore, *positionAfter;
 int nPlanes;

} TRANSECTMORPH_PARMS;

void * transectMorph_parms(TRANSFORM * t);
void   transectMorph(TRANSFORM * t);

#endif 
