#include "recordDecorator.h"
#include <string.h>
#include "ddcMalloc.h"

RECORD_DECORATOR* recordDecorator_init(
   unsigned lrec, unsigned nfields,
   const char* fieldTypes, const char* fieldNames, const char* fieldUnits,
   rd_printFunction_t printFunction)
{
   RECORD_DECORATOR* this = ddcMalloc(sizeof(RECORD_DECORATOR));

   this->lrec = lrec;
   this->nfields = nfields;
   this->fieldTypes = strdup(fieldTypes);
   this->fieldNames = strdup(fieldNames);
   this->fieldUnits = strdup(fieldUnits);
   this->print = printFunction;
   return this;
}
