#ifndef BACK_COMMUNICATE__
#define BACK_COMMUNICATE__

#include <mpi.h>
#include "integer.h"
#include "ddc.h"

typedef struct {
  int left,right,parent,rank;
  double c0[3],c1[3];
} nodedata;

typedef struct domain_info_s *domain_info_p;

domain_info_p redistribute2(integer np,integer pid,MPI_Comm comm,
			    integer nloc,PARTICLE pvec[],
			    double c0[3],double c1[3]);


integer dinf_get_neilist(domain_info_p dip,double rcut,nodedata **list);
const nodedata * dinf_get_self(domain_info_p dip);
const nodedata * dinf_lookup_point(domain_info_p dip,double x[3]);

void dinf_printtree(domain_info_p dip);
void dinf_delete(domain_info_p dip);

int dinf_overlap(double a0[3],double a1[3],double b0[3],double b1[3]);
int dinf_inside(const double c0[3],const double c1[3],const double x[3]);

integer treesize(integer n);
integer treeself(integer p0,integer p1,integer p);


#endif
