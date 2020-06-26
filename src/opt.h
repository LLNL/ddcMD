#ifndef OPT_H
#define OPT_H

#include "gid.h"
#include "three_algebra.h"

#define OPT_PARTICLE THREE_VECTOR 

typedef struct opt_nbox_st { int offset , image; } OPT_NBOX;
typedef struct opt_box_st { int first, last, nlocal, nremote, nn, ni, edge, image ; int x,y,z; OPT_NBOX *nlist; OPT_NBOX *ilist; } OPT_BOX; 
typedef struct opt_box_size_st { int nx,ny,nz,nbox,nimage; } OPT_BOX_SIZE; 

typedef struct opt_particle0_st { double x,y,z; int type; gid_type label;} OPT_PARTICLE0; 

typedef struct opt_strip_st {int last,first,image;} OPT_STRIP ; 

// not implemented int opt_sort_by_box(OPT_PARTICLE *a, OPT_PARTICLE *b);
// not implemented int opt_sort_by_image(OPT_STRIP *a, OPT_STRIP *b);
OPT_BOX *opt_box(double rcut, int n, int nlocal, OPT_PARTICLE *rv, unsigned *sort_index, OPT_BOX_SIZE *box_size);
int opt_transform(double rcut,int n, int nlocal, double *rx, double *ry, double *rz, OPT_PARTICLE *rv) ;
int opt_define_strips(OPT_BOX *box, OPT_STRIP *strip) ;
int opt_countNeighbors(OPT_BOX *box);

#endif // #ifndef OPT_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
