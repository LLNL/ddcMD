#include "collection.h"
#include "state.h"

#ifndef MINMAXGPU_H_
#define MINMAXGPU_H_
#define SHARED_RED 64
#ifdef __cplusplus
extern "C"
{
#endif
double* getMinMaxGPU(COLLECTION *col, STATE * gpuState_h, STATE* hostState, double *mmgpu);
void execReduceXYZ(double *g_idata_x, double *g_idata_y, double *g_idata_z,
                   double *g_odata_xMax, double *g_odata_yMax, double *g_odata_zMax,
                   double *g_odata_xMin, double *g_odata_yMin, double *g_odata_zMin, double *mmgpu, int n);
double *getMax(double *in, int n);
long timeMinMaxThrust(int len);

#ifdef __cplusplus
}
#endif

#endif //MINMAXGPU_H_
