#ifndef RECORD_DECORATOR_H
#define RECORD_DECORATOR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*rd_printFunction_t)(char* buf, int recordId);

/** Supports functionality that allows objects to register additional
 *  fields in records for output files with RecordDecorator support.  If
 *  this were C++ then the print function would be a functor.  As it is,
 *  clients like voronoiLoadBalance end up using static local variables
 *  to store print state.  There are other solutions available, even in
 *  C, that would allow us to avoid static variables, but they incur
 *  other complications.  For now, we'll stick with simplicity and
 *  develop something more elegant if it really becomes necessary.
 */
typedef struct RecordDecorator_st
{
   unsigned lrec;
   unsigned nfields;
   char* fieldTypes;
   char* fieldNames;
   char* fieldUnits;

   rd_printFunction_t print;
} RECORD_DECORATOR;

RECORD_DECORATOR* recordDecorator_init(
   unsigned lrec, unsigned nfields,
   const char* fieldTypes, const char* fieldNames, const char* fieldUnits,
   rd_printFunction_t printFunction);

#ifdef __cplusplus
}
#endif

#endif
