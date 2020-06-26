#include "collection.h"
#include "three_algebra.h"
#include "expandbuffer.h"
#include "error.h"

typedef struct trvf_struct { unsigned type; THREE_VECTOR r,v,f;} TRVF; 
static THREE_VECTOR *__v = NULL; 
static TRVF *__trvf = NULL; 
THREE_VECTOR *state_getVPtr(STATE *state)
{
	__v = ExpandBuffers((void *)__v, sizeof(THREE_VECTOR), state->nion, 512, LOCATION("state_getVPtr"), "__v");
	for (int i=0;i<state->nion;i++) 
	{
		__v[i].x = state->vx[i];
		__v[i].y = state->vy[i];
		__v[i].z = state->vz[i];
	}
	return __v; 
}
THREE_VECTOR *state_getTrvfPtr(STATE *state)
{
	__v = ExpandBuffers((void *)__trvf, sizeof(TRVF), state->nion, 512, LOCATION("state_getTrvfPtr"), "__Trvf");
	for (int i=0;i<state->nion;i++) 
	{
		__trvf[i].type=state->atomtype[i] & 0xffff  ;
		__trvf[i].r.x = state->rx[i];
		__trvf[i].r.y = state->ry[i];
		__trvf[i].r.z = state->rz[i];
		__trvf[i].v.x = state->vx[i];
		__trvf[i].v.y = state->vy[i];
		__trvf[i].v.z = state->vz[i];
		__trvf[i].f.x = state->fx[i];
		__trvf[i].f.y = state->fy[i];
		__trvf[i].f.z = state->fz[i];
	}
	return __v; 
}
void state_putV(STATE *state,THREE_VECTOR __v,int i)
{
		state->vx[i]=__v.x ;
		state->vy[i]=__v.y ;
		state->vz[i]=__v.z ;
}
void state_putTrvfPtr(STATE *state,TRVF trvf,int i)
{
	for (int i=0;i<state->nion;i++) 
	{
		state->atomtype[i] = trvf.type+(state->atomtype[i] & 0xffff0000)  ;
		state->rx[i]=trvf.r.x;
		state->ry[i]=trvf.r.y;
		state->rz[i]=trvf.r.z;
		state->vx[i]=trvf.v.x;
		state->vy[i]=trvf.v.y;
		state->vz[i]=trvf.v.z;
		state->fx[i]=trvf.f.x;
		state->fy[i]=trvf.f.y;
		state->fz[i]=trvf.f.z;
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
