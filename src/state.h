#ifndef STATE_H
#define STATE_H
#include "three_algebra.h"
#include "species.h"
#include "group.h"
//#include "random.h"
typedef struct state_st
{
	double *rx, *ry, *rz; 
   double *vx, *vy, *vz; 
   double *fx, *fy, *fz; 
   double *q; 
   double *potentialEnergy;
   double *virial;
	THREE_SMATRIX *sion; 
	gid_type* label;
	int *atomtype;
	SPECIES **species;
	struct group_st  **group;
	int nlocal;
	int nion;

  float *cost;

  int domcen_active,domcen_np,domcen_pid;
  double (*domcen_vec)[3];
} STATE;
THREE_VECTOR *state_getVPtr(STATE *state);
void state_putV(STATE *state,THREE_VECTOR __v,int i);
void resize(int n, int mode, STATE*state);
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
