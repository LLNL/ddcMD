#ifndef MD2DDC_H
#define MD2DDC_H
#include <mpi.h>
#include "gid.h"

void *GetLocationinfo(int nRecv, int RecvStart, MPI_Datatype *RECVTYPE);
void *GetRecvLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE);
void *GetVelocityLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE);
void *GetForceLocation(int n, int Start, MPI_Datatype*TYPE,unsigned calculate);
void fillsendbuf(double *send, int nlist, int *list);
void Velocityfillsendbuf(double *send, int nlist, int *list);
void forceaccum0(double *force, int nlist, int *list);
void forceaccum1(double *force, int nlist, int *list);
void forceaccum2(double *force, int nlist, int *list);
void forceaccum3(double *force, int nlist, int *list);
void forceaccum4(double *force, int nlist, int *list);
void forceaccum5(double *force, int nlist, int *list);
void forceaccum6(double *force, int nlist, int *list);
void forceaccum7(double *force, int nlist, int *list);
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
