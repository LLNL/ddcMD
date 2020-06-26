#ifndef PROCESSGRID2D_H
#define PROCESSGRID2D_H

#include <mpi.h>

typedef struct pgrid_data_struct 
{
  MPI_Comm comm_;
  int nprow_,npcol_;
} PGRID_DATA;

void pgrid_init(PGRID_DATA *pdata, MPI_Comm comm, int nprow,int npcol);
int pgrid_mype(PGRID_DATA *pdata);
int pgrid_npes(PGRID_DATA *pdata);
int pgrid_myrow(PGRID_DATA *pdata, int mype);
int pgrid_mycol(PGRID_DATA *pdata, int mype);
int pgrid_pe(PGRID_DATA *pdata, int i, int j);
MPI_Comm pgrid_comm(PGRID_DATA *pdata);
MPI_Comm pgrid_subcomm(PGRID_DATA *pdata, int nrow, int ncol,int istart, int jstart);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
