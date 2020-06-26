#include "saveState.h"
#include "gid.h"
#include "error.h"     // LOCATION
#include "ptiming.h"     // LOCATION
#include "io.h"
#include "expandbuffer.h"

#include "mpiUtils.h"
#include <assert.h>
#include "three_algebra.h"

THREE_MATRIX box_get_h(BOX_STRUCT *box);          // defined in box.c
void box_put(BOX_STRUCT*box, int cmd, void *ptr); // defined in box.c

typedef struct stateBackup_st
{
   double _time;
   int64_t _loop;
   int _nlocal;
   SIMULATE* _simulate;
   THREE_MATRIX _h;
   gid_type* _label;
   int* _atype; // atomtype
   SPECIES** _species;
   GROUP** _group;
   double* _rx;
   double* _ry;
   double* _rz;
   double* _vx;
   double* _vy;
   double* _vz;
   DOMAINX* _domainx;
} StateBackup;

static StateBackup _b = {0.0,0,0,NULL,
                         {0.0,0.0,0.0,
                          0.0,0.0,0.0,
                          0.0,0.0,0.0},
                         NULL,NULL,NULL,NULL,
                         NULL,NULL,NULL,NULL,NULL,NULL,
                         NULL} /* Other elements (mem)set to 0 by C language spec. */;

/** This function may be called by tasks that have simulate==NULL (i.e.,
 *  slave tasks) */
int saveState(SIMULATE* simulate)
{
/*
  if(getRank(0) == 0)
    printf("In %s() at %s:%d: State saving not really supported any more, "
           "but you can uncomment this abort to try:-)\n",
           __func__,__FILE__,__LINE__);
  abortAll(13);
*/

   if (simulate == NULL) return 0;
   
   profile(SAVESTATE, START);
   SYSTEM* sys = simulate->system;
   STATE* state = simulate->system->collection->state;
   DOMAINSET* domainset = &(simulate->ddc->domains);

   /*
     The 425000 was for BG/L time, and buffer set as big as
      needed from start to avoid memory fragmentation. 5000
      seems more reasonable.
   */
   int incr=5000 /*425000*/ ;

   int nlocal = sys->nlocal;
   const int dblsz = sizeof(double);
   const int intsz = sizeof(unsigned);
   const int ptrsz = sizeof(void*);
   const int dmnxsz = sizeof(DOMAINX);
   
   _b._label = ExpandBuffers((void*) _b._label, sizeof(gid_type),nlocal, incr, LOCATION("saveState"), "label_backup");
   _b._atype = ExpandBuffers((void*) _b._atype, intsz, nlocal, incr, LOCATION("saveState"), "atomtype_backup");
   _b._species = ExpandBuffers((void*) _b._species, ptrsz, nlocal, incr, LOCATION("saveState"), "species_backup");
   _b._group = ExpandBuffers((void*) _b._group, ptrsz, nlocal, incr, LOCATION("saveState"), "group_backup");
   _b._rx = ExpandBuffers((void*) _b._rx, dblsz, nlocal, incr, LOCATION("saveState"), "rx_backup");
   _b._ry = ExpandBuffers((void*) _b._ry, dblsz, nlocal, incr, LOCATION("saveState"), "ry_backup");
   _b._rz = ExpandBuffers((void*) _b._rz, dblsz, nlocal, incr, LOCATION("saveState"), "rz_backup");
   _b._vx = ExpandBuffers((void*) _b._vx, dblsz, nlocal, incr, LOCATION("saveState"), "vx_backup");
   _b._vy = ExpandBuffers((void*) _b._vy, dblsz, nlocal, incr, LOCATION("saveState"), "vy_backup");
   _b._vz = ExpandBuffers((void*) _b._vz, dblsz, nlocal, incr, LOCATION("saveState"), "vz_backup");
   _b._domainx = ExpandBuffers((void*) _b._domainx, dmnxsz, domainset->size, 16, LOCATION("saveState"), "domainx_backup");
   
   _b._time     = sys->time;
   _b._loop     = sys->loop;
   _b._nlocal   = sys->nlocal;
   _b._simulate = simulate;
   _b._h        = box_get_h(sys->box);

   for (int ii=0; ii<nlocal; ++ii)
   {
      _b._label[ii] = state->label[ii];
      _b._atype[ii] = state->atomtype[ii];
      _b._species[ii] = state->species[ii];
      _b._group[ii] = state->group[ii];
      _b._rx[ii] = state->rx[ii];
      _b._ry[ii] = state->ry[ii];
      _b._rz[ii] = state->rz[ii];
      _b._vx[ii] = state->vx[ii];
      _b._vy[ii] = state->vy[ii];
      _b._vz[ii] = state->vz[ii];
   }

   // save domain centers and radii
   for (int ii=0; ii<domainset->size; ++ii)
      _b._domainx[ii] = domainset->domains[ii];
	
   profile(SAVESTATE, END);
   return 0;
}

/** This function may be called by tasks that have simulate==NULL (i.e.,
 *  slave tasks) */
int restoreState(SIMULATE* simulate)
{

   if (! simulate) return 0;

   SYSTEM* sys = simulate->system;
   STATE* state = simulate->system->collection->state;
   int nlocal = _b._nlocal;
   
   sys->time = _b._time;
   sys->loop = _b._loop;
   sys->nlocal = _b._nlocal;
   box_put(sys->box, HO, &(_b._h));
   box_put(sys->box, HFAC,  (void *)&I_3x3);
   
   simulate->time = _b._time;
   simulate->loop = _b._loop;
      
   for (int ii=0; ii<nlocal; ++ii)
   {
      state->label[ii]    = _b._label[ii];
      state->atomtype[ii] = _b._atype[ii];
      state->species[ii]  = _b._species[ii];
      state->group[ii]    = _b._group[ii];
      state->rx[ii]       = _b._rx[ii];
      state->ry[ii]       = _b._ry[ii];
      state->rz[ii]       = _b._rz[ii];
      state->vx[ii]       = _b._vx[ii];
      state->vy[ii]       = _b._vy[ii];
      state->vz[ii]       = _b._vz[ii];
   }

     DOMAINSET* domainset = &(simulate->ddc->domains);
     for (int ii=0; ii<domainset->size; ++ii)
       domainset->domains[ii] = _b._domainx[ii];
   return 0;
}

int
writeState(void)
{
   restoreState(_b._simulate);

   CreateSnapshotdir(_b._simulate, NULL);
   writeRestart(_b._simulate,1);

   return 0;
}


/* Local Variables: */
/* tab-width: 3 */
/* indent-tabs-mode: nil */
/* End: */
