#include "qIntrinsics.h"

extern inline double * force_align_4(double *p);
extern inline double * force_align_2(double *p);

#ifndef __VECTOR4DOUBLE__
extern inline v4d vec_lds(long off,double *p);
extern inline v4d vec_ld(long off,double *p);
extern inline void vec_st(v4d v,long off,double *p);
extern inline void vec_st2(v4d v,long off,double *p);
extern inline v4d vec_ld2(long off,double *p);
extern inline v4d vec_splats(double x);
extern inline v4d vec_splat(v4d x,int k);
extern inline v4d vec_mul(v4d x,v4d y);
extern inline v4d vec_add(v4d x,v4d y);
extern inline v4d vec_sub(v4d x,v4d y);
extern inline v4d vec_madd(v4d x,v4d y,v4d z);
extern inline v4d vec_msub(v4d x,v4d y,v4d z);
extern inline v4d vec_nmsub(v4d x,v4d y,v4d z);
extern inline v4d vec_sldw(v4d x,v4d y,int k);
extern inline v4d vec_gpci(int x);
extern inline v4d vec_perm(v4d a,v4d b,v4d p);
extern inline v4d vec_lvsl(long off,double *p);
extern inline v4d vec_cfidu(v4d x);
extern inline v4d vec_swdiv(v4d num,v4d denom);
extern inline v4d vec_swdiv_nochk(v4d num,v4d denom);
extern inline v4d vec_re(v4d x);
extern inline double vec_extract(v4d x,int k);
#endif
