#ifndef NPTGLF_H
#define NPTGLF_H
typedef struct nptglf_parms_st
{
	double tau, Gamma, zeta, pressure;
	int randomState[35];
} NPTGLF_PARMS;
#endif
