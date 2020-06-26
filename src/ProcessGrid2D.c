//
// ProcessGrid2D:  functions to manage a 2D process grid of nprow x npcol processors.
//
// written by Erik Draeger, LLNL, 12/19/2008


#include <mpi.h>
#include "ProcessGrid2D.h"

////////////////////////////////////////////////////////////////////////////////
void pgrid_init(PGRID_DATA *pdata, MPI_Comm comm, int nprow, int npcol) {
  pdata->comm_ = comm;
  pdata->nprow_ = nprow;
  pdata->npcol_ = npcol;
}
////////////////////////////////////////////////////////////////////////////////
int pgrid_mype(PGRID_DATA *pdata) {
  int mype;
  MPI_Comm_rank(pdata->comm_, &mype);
  return mype;
}
////////////////////////////////////////////////////////////////////////////////
int pgrid_npes(PGRID_DATA *pdata) {
  int npes;
  MPI_Comm_size(pdata->comm_, &npes);
  return npes;
}
////////////////////////////////////////////////////////////////////////////////
int pgrid_myrow(PGRID_DATA *pdata, int pe) {
  int px = pe;
  while (px >= pdata->nprow_) px -= pdata->nprow_;
  return px;
}
////////////////////////////////////////////////////////////////////////////////
int pgrid_mycol(PGRID_DATA *pdata, int pe) {
  int py = pe/pdata->nprow_;
  return py;
}
////////////////////////////////////////////////////////////////////////////////
int pgrid_pe(PGRID_DATA *pdata, int i, int j) {
  int ptmp = i + j*pdata->nprow_;
  return ptmp;
}
////////////////////////////////////////////////////////////////////////////////
MPI_Comm pgrid_comm(PGRID_DATA *pdata) {
  return pdata->comm_;
}
////////////////////////////////////////////////////////////////////////////////
MPI_Comm pgrid_subcomm(PGRID_DATA *pdata, int npex, int npey, int istart, int jstart) {
  //assert(npex <= nprow_);
  //assert(npey <= npcol_);

  int nprocs = npex*npey;
  int pmap_[nprocs];
  // build pmap
  int p = 0;

  int i,j;
  for (j=jstart; j<npey+jstart; j++) {
    for (i=istart; i<npex+istart; i++) {
      pmap_[p] = i + j*pdata->nprow_;
      p++;
    }
  }

  MPI_Comm subcomm_;
  MPI_Group c_group, subgroup;
  MPI_Comm_group(pdata->comm_,&c_group);
  MPI_Group_incl(c_group,nprocs,&pmap_[0],&subgroup);
  MPI_Comm_create(pdata->comm_,subgroup,&subcomm_);
  MPI_Group_free(&c_group);
  MPI_Group_free(&subgroup);
  return subcomm_;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
