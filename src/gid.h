#ifndef GID_H
#define GID_H
#include <inttypes.h>
#define __STDC_FORMAT_MACROS
//#include <stdint.h>
#define LONG64  uint64_t
#define SIGNED64  int64_t
//#define LONG64 long long unsigned
#define MPI_GID_TYPE MPI_UNSIGNED_LONG_LONG
typedef LONG64 gid_type;
static const gid_type gid_max=(~((gid_type )0))-1;
static const gid_type invalidGid = 0xffffffffffffffff; 
//static const char * invalidGidString = "0xffffffffffffffff"; 
#define invalidGidString  "0xffffffffffffffff" 
int compareGid(const void *v1, const void *v2);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
