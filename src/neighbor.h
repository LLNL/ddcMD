#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "gid.h"   // LONG64
#include "geom.h"  // PARTICLESET

#define MAXNC    4
#define NO_NBR_UPDATE   -174321321


/**
 *  The neighbor list has the following responsibilities:
 *
 *  - Provide a (one sided) list of all pairs that can be iterated from
 *  multiple cutoffs.
 *
 *  - Regenerate the pair list on demand
 *
 *  - Detect when the pair list needs to be rebuilt due to either motion
 *  of the particles or re-assignment of particles to tasks.
 *
 *  - Update stored properties of pairs that change every time step
 *  (such as pair distances).
 */



/** These modes determine what kinds of bonds will be kept in the nbr
 * list for a given cutoff.
 *
 * - RCUT_NONE ->   Not used.
 * - RCUT_LOCAL ->  First vertex of any pair will be a local particle.
 * - RCUT_REMOTE -> Not used.
 * - RCUT_ALL ->    First vertex of any pair can be either local or
 *                  remote atom.
 *
 * In practice, as of r630, if rcut[0].mode == RCUT_LOCAL the nbr code
 * will exclude remote-local and remote-remote pairs from the list (only
 * local-local and local-remote pairs will be in the list).  Any other
 * setting will allow all four types of pairs.
 */ 
enum RCUT_ENUMS { RCUT_NONE=0, RCUT_LOCAL=1, RCUT_REMOTE=2, RCUT_ALL=3 } ;
enum NEIGHBORTABLETYPE { NEIGHBORTABLE_NONE=0, NEIGHBORTABLE_SKINNY=1, NEIGHBORTABLE_FAT=2,  NEIGHBORTABLE_GPU=4};

typedef struct rcut_str
{
   double value;
   enum RCUT_ENUMS  mode;
   int type;
} RCUT_TYPE; 

/** hmat was designed for mgpt, but can be thought of as potential
 *  dependent pairwise storage.  Any potential can just ensure that each
 *  rpair points at an appropriate chunk of memory. */
typedef struct rpair_struc
{
   double x, y, z, r;
   THREE_VECTOR fp;   // accumulation for virtual bond forces (VBF).
   double e;          // accumulation for particle energy when using VBF.
   void *hmat;
} RPAIR;

typedef struct p2list_struc
{
	RPAIR *q;
	int j, flag;
} P2LIST;

typedef struct p3list_struc
{
	RPAIR *q[3];
	int klocal, flag;
} P3LIST;

typedef struct pairs_struc
{
/*	long long unsigned int pkey; */
	int j; 
	struct pairs_struc *ilink;
	RPAIR *p;
} PAIRS;

typedef struct nparticle_struc
{
	THREE_VECTOR r0;
	int nb;
	PAIRS *ifirst[MAXNC], *ilast[MAXNC];
	int tmp; 
} NPARTICLE;


typedef struct tlist_struc
{
	int ij, jk, ik, i, j, k;
	int link, r;
} TLIST;

typedef struct qlist_struc
{
	int ij, jk, kl, li, ik, jl, flag;
} QLIST;

/**
 * update is set but the value is never used.
 */
typedef struct nbr_struc
{
	char *name;		/* name of the system */
	char *objclass;
	char *value;
	char *type;		/* type label */
	void *parent; 
	unsigned nSearch, npairs, nlocal, number_particles;
	unsigned nRpair, update, pair_per_particle_hint;
	int nc; // number of cutoffs
	THREE_VECTOR center,rbar;
	NPARTICLE *particles;
	THREE_MATRIX h0; 
	THREE_MATRIX h0_inv; 
	double rcut0;
	double deltaR;
	RCUT_TYPE* rcut;
	GEOM *geom;
	PAIRS *pairs;
	RPAIR *Rpairs;
   unsigned Rpairs_blk;
   int lastUpdate;
} NBR;

int neighborCheck(NBR *nbr, PARTICLESET*particleset);
void neighborCutoff(NBR *nbr,int n, RCUT_TYPE *rcut);
void neighborRef(NBR *nbr, PARTICLESET*particleset);
NBR *neighbors(NBR *nbr, PARTICLESET*particleset);
NBR *neighbors1(NBR *nbr, PARTICLESET*particleset);
void pairUpdate(NBR *nbr, PARTICLESET*particleset);
void pairelement(NBR* nbr, NPARTICLE *pi, PAIRS *pij, RPAIR *p, double x, double y, double z,double r2,int local_remote_flag);
#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
