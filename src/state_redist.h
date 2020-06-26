#ifndef STATE_REDIST__
#define STATE_REDIST__

#include <mpi.h>
#include "integer.h"
#include "state.h"

integer state_redist(STATE *state_p,MPI_Comm comm,integer nitems);

#endif
