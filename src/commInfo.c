#include "commInfo.h"

#ifdef BGP
#include <mpix.h>
#include <mpido_properties.h>
int isRectangularCommBgp(MPI_Comm comm);
#endif



/** Returns non-zero for rectangular comms. Returns zero for
 * MPI_COMM_NULL.  Returns zero on machines with no notion of what a
 * rectangular comm is. */
int isRectangularComm(MPI_Comm comm)
{
   if (comm == MPI_COMM_NULL)
      return 0;
   
   
   #ifdef BGP
   return isRectangularCommBgp(comm);
   #endif

   return 0;
}

#ifdef BGP
int isRectangularCommBgp(MPI_Comm comm)
{
   int value;
   MPIX_Get_property(comm, MPIDO_RECT_COMM, &value);
   return value;
}
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
