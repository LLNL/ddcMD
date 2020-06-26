#ifndef PREFIX_SCAN_H_
#define PREFIX_SCAN_H_
#define SHARED_PREFIX_SCAN 128

#ifdef __cplusplus
extern "C"{ 
#endif

void plusScan(double *in, double *out, double *scratch, double *scratch1, int n );
void plusScanI(int *in, int *out, int *scratch, int *scratch1, int n );

#ifdef __cplusplus
}

#endif
#endif
