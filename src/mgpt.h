#ifndef MGPT_H 
#define MGPT_H 
#include "neighbor.h"
#include "mgpt5X5.h"
#include "mgpt7X7.h"
typedef struct hmatinfo_st
{
	int order, size, sizeof_hmat;
	double s, p, d, f, anorm3, anorm4;
	int nps;
	void *hmat;
	double (*AB) ();
	double (*V33FM) ();
	double (*V4) ();
	void (*HAMLTN) ();
} HAMLTNINFO;
typedef struct v2parm_str
{
	double vpair[4], ktan[4], dvdvol[4];
}
V2PARM;

typedef struct mgptparms_st
{
	double vol, rmax, rpair, rcrit;
	double va, vb, vc, vd, ve, pa, pb, pc, pd, pe, evol0, pvol0;
	double deltaR,recipDeltaR; 
	double *v2parm; 
	double *rmesh; 
	double *xKnots; 
	double *volKnots; 
	double vol0; 
	double r0rws,rcrws,rmrws; 
	double p1,dp1; 
	double al,dal; 
	double alm,dalm; 
	int mode,ipot; 
	int nRMesh, nVolumes;
	int nvmin,nvmax; 
	int volume_from_box;
    int debug;
	LONG64  TempStorageSize; 
	HAMLTNINFO *hamltninfo;
	char *parmfile; 
} MGPT_PARMS;

HMATD *hamltnd(RPAIR*r, int mode);
HMATF *hamltnf(RPAIR*r, int mode);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
