#ifndef MD_H 
#define MD_H 
typedef struct mdparm_st
{
	int istep, iprint, ilist, ascii_restart, binary_restart, ascii_xyz, volume_control;
	double dt, tau, Gamma, csi, c2;
	double pstart, pend, pt0, prate;
	double tstart, tend, tt0, trate;
} MDPARM;
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
