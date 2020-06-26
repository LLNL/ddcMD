#ifndef RECNLIST__
#define RECNLIST__

#include <mpi.h>
#include "integer.h"

typedef struct {
  double xx0[3],xx1[3];
  int pid0,pid1,stride,padding;
} recnlist_node;


void recnlist_refine(integer pid,MPI_Comm comm,
		     recnlist_node *self,recnlist_node *sibling,double rcut,
		     int nnei_p[1],recnlist_node **nlist_p);

#endif
