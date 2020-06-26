#ifndef BGQ_INTRINSICS__
#define BGQ_INTRINSICS__

//#ifndef BGQ
#ifdef __VECTOR4DOUBLE__
/* If we have the builtin vector intrinsics (e.g. XL compiler on BG/Q) */
#include <stdlib.h>
typedef size_t ptrint_t;

typedef vector4double v4d;

inline double * force_align_4(double *p) {
  ptrint_t i = (ptrint_t) p;
  return (double *) (i - (i&31));
}

inline double * force_align_2(double *p) {
  ptrint_t i = (ptrint_t) p;
  return (double *) (i - (i&15));
}

#else
/* Implmenetation of many BG/Q vector intrinsics */
#include <assert.h>

#include <stdlib.h>
typedef size_t ptrint_t;

typedef struct {
  double v[4];
} vector4double;
typedef vector4double v4d;

inline v4d vec_lds(long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = poff[0];
  return v;
}

inline double * force_align_4(double *p) {
  ptrint_t i = (ptrint_t) p;
  //assert((i&31) == 0);
  return p;
  return (double *) (i - (i&31));
}

inline double * force_align_2(double *p) {
  ptrint_t i = (ptrint_t) p;
  //assert((i&15) == 0);
  return p;
  return (double *) (i - (i&15));
}

inline v4d vec_ld(long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);
  v4d v;
  double *q = force_align_4(poff);
  int i;
  for(i = 0; i<4; i++) v.v[i] = q[i];
  return v;
}
inline void vec_st(v4d v,long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);  
  double *q = force_align_4(poff);
  int i;
  for(i = 0; i<4; i++) q[i] = v.v[i];
}
inline void vec_st2(v4d v,long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);  
  double *q = force_align_2(poff);
  int i;
  for(i = 0; i<2; i++) q[i] = v.v[i];
}

inline v4d vec_ld2(long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);
  double *q = force_align_2(poff);
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = q[i&1];
  return v;
}

inline v4d vec_splats(double x) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x;
  return v;
}

inline v4d vec_splat(v4d x,int k) {
  v4d v;
  int i;
  assert(k >= 0 && k < 4);
  for(i = 0; i<4; i++) v.v[i] = x.v[k];
  return v;
}

inline v4d vec_mul(v4d x,v4d y) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x.v[i]*y.v[i];
  return v;
}

inline v4d vec_add(v4d x,v4d y) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x.v[i] + y.v[i];
  return v;
}

inline v4d vec_sub(v4d x,v4d y) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x.v[i] - y.v[i];
  return v;
}
inline v4d vec_madd(v4d x,v4d y,v4d z) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x.v[i]*y.v[i] + z.v[i];
  return v;
}
inline v4d vec_msub(v4d x,v4d y,v4d z) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = x.v[i]*y.v[i] - z.v[i];
  return v;
}
inline v4d vec_nmsub(v4d x,v4d y,v4d z) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = -(x.v[i]*y.v[i] - z.v[i]);
  return v;
}

inline v4d vec_sldw(v4d x,v4d y,int k) {
  v4d v;
  int i;
  assert(k >= 0 && k < 4);
  for(i = 0; i<4; i++)
    v.v[i] = (k+i < 4) ? x.v[k+i] : y.v[k+i-4];
  return v;
}

inline v4d vec_gpci(int x) {
  v4d v;
  int i;
  assert(x >= 0 && x < (1<<(4*3)));
  for(i = 0; i<4; i++)
    v.v[i] = 2.0 + 0.25*((x>>(3*(3-i))) & 7);
  return v;
}
inline v4d vec_perm(v4d a,v4d b,v4d p) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) {
    int idx = (int) ((p.v[i]-2.0)*4.0);
    assert(p.v[i] == (idx*0.25+2.0) && idx>=0 && idx<8);
    v.v[i] = (idx < 4) ? a.v[idx] : b.v[idx-4];
  }
  return v;
}
inline v4d vec_lvsl(long off,double *p) {
  double *poff = (double *) (((ptrint_t) p) + off);  
  double *palign = force_align_4(poff);
  const int d = poff - palign;
  int i,perm = 0;
  assert(d >= 0 && d<4);
  for(i = 0; i<4; i++)
    perm = (perm<<3) | (d+i);
  return vec_gpci(perm);
}
inline v4d vec_cfidu(v4d x) {
  int i;
  union {
    v4d vec;
    unsigned long long int i64[4];
  } conv;
  assert(sizeof(double) == sizeof(unsigned long long int) &&
	 sizeof(conv) == sizeof(v4d));
  conv.vec = x;
  for(i = 0; i<4; i++)
    conv.vec.v[i] = (double) conv.i64[i];
  return conv.vec;
}
inline v4d vec_swdiv(v4d num,v4d denom) {
  v4d v;
  int i;
  for(i = 0; i<4; i++)
    v.v[i] = num.v[i] / denom.v[i];
  return v;
}
inline v4d vec_swdiv_nochk(v4d num,v4d denom) {
  return vec_swdiv(num,denom);
}
inline v4d vec_re(v4d x) {
  v4d v;
  int i;
  for(i = 0; i<4; i++) v.v[i] = 1.0/x.v[i];
  return v;
}
inline double vec_extract(v4d x,int k) {
  assert(k >= 0 && k < 4);
  return x.v[k];
}

#endif

#endif
