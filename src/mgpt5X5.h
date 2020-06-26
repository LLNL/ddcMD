#ifndef MGPT5X5_H 
#define MGPT5X5_H 
#include "gid.h"
#define SIZE_M5X5  25
#define SIZE_C5    5
// The matrixes v,x,y,z are double word aligned 

typedef struct hmatd_st
{
	double v[25], f, x[25], df,y[25], df_dvol,z[25];
	int cnt;
	unsigned short min, max;
} HMATD;
#if defined(BGL) || defined(BGP) 

#define AB5X5             AB5X5_bg
#define TABC5X5X3         TABC5X5X3_bg
#define TABC5X5X4         TABC5X5X4_bg

#endif 

#if !defined(AB5X5) 

#define AB5X5  AB5X5_1
#define TABC5X5X3  TABC5X5X3_1
#define TABC5X5X4  TABC5X5X4_1
#endif 

#endif 



/* Local Variables: */
/* tab-width: 3 */
/* End: */
