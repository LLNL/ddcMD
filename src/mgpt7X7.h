#ifndef MGPT7X7_H 
#define MGPT7X7_H 
//#define HMATF_PADDED
#define HMATF_TWISTED

#ifdef HMATF_PADDED
    typedef struct hmatf_st
    {
	   double v[56], x[56], y[56], z[56],f,df,df_dvol;
	   int cnt;
	   int scnt;
	   //unsigned short min, max;
    } HMATF;

#   define SIZE_M7X7  56
#   define SIZE_C7    8
#   define HMATF_FORMAT padded
#if defined(BGL) || defined(BGP)
#     define AB7X7  AB7X7_bg_padded
#	  define TABC7X7X3  TABC7X7X3_bg_padded
#	  define TABC7X7X4  TABC7X7X4_bg_padded
#   endif 
#   ifndef AB7X7 
# 	  define AB7X7  AB7X7_padded
#	  define TABC7X7X3  TABC7X7X3_padded
#	  define TABC7X7X4  TABC7X7X4_padded
#   endif 
#	define HMATF_STORE  storeRegular
#endif 

#ifdef HMATF_TWISTED
#   define SIZE_M7X7  56
#   define SIZE_C7    8
#   define HMATF_FORMAT twisted
   typedef struct hmatf_st
   {
	   double v[56], f, x[56], df, y[56], df_dvol, z[56];
	   int cnt;
	   unsigned short min, max;
   } HMATF;
#if defined(BGL) || defined(BGP)
#	  ifndef DEVLPMT
#        define AB7X7  AB7X7_bg_twisted
#	     define TABC7X7X3  TABC7X7X3_bg_twisted
#	     define TABC7X7X4  TABC7X7X4_bg_twisted
#     else
#        define AB7X7  AB7X7_bg_twisted_devlpmt
#	     define TABC7X7X3  TABC7X7X3_bg_twisted_devlpmt
#	     define TABC7X7X4  TABC7X7X4_bg_twisted_devlpmt
#     endif 
	
#endif 
#ifndef AB7X7 
#	define AB7X7  AB7X7_twisted
#	define TABC7X7X3  TABC7X7X3_twisted
#	define TABC7X7X4  TABC7X7X4_twisted
#endif 
#	define HMATF_STORE  storeTwisted
#endif 

#ifndef HMATF_FORMAT 
#   define SIZE_M7X7  49
#   define SIZE_C7    7
   typedef struct hmatf_st
   {
	   double v[49], f, x[49], df, y[49], df_dvol, z[49];
	   int cnt;
	   unsigned short min, max;
   } HMATF;

#if defined(BGL) || defined(BGP)
#     define AB7X7  AB7X7_bg
#	  define TABC7X7X3  TABC7X7X3_bg
#	  define TABC7X7X4  TABC7X7X4_bg
#endif 
#ifndef AB7X7 
#     define AB7X7  AB7X7_1
#	  define TABC7X7X3  TABC7X7X3_1
#	  define TABC7X7X4  TABC7X7X4_1
#endif
#define HMATF_STORE  storeRegular
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
