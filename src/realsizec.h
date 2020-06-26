#ifndef REALSIZEC_H
#define REALSIZEC_H
#include "codata.h"
/*
  don't use John's atomic scale units:
  time:  1.460509 x 10^(-15) sec
  distance: Bohr [0.52917 x 10^(-10) m]
  energy: Rydbergs [ 2.1798741 x 10^(-18) J = 13.605698 eV]
  mass: amu [1.66054 x  10^(-27) kg]
  pressure: 14711.12 GPa = 145.1875 Mbar ]

  instead, use time in fs and mass in 0.4688u   0.468791526u
*/
#define USE_CODATA
#ifdef USE_CODATA
#define amu2M ( u_MKS * 1e+30*(a0_MKS*a0_MKS)/Rinfhc_MKS ) 
#define Ry2eV Rinfhc_eV         
#define au2a ( a0_MKS * 1e10 )
#define Ry2K  (Rinfhc_eV / kB_eV)
#define P2GPa Rinfhc_MKS/(a0_MKS*a0_MKS*a0_MKS) * 1e-9 
#define re  (re_MKS* 1e10*a2au ) 
#define mec2  ((mec2_MeV * 1e6)/ Eh_eV )
#define hc   ( 1e10 *a2au/Rinf )
#else 
#define amu2M 2.1331426 
#define Ry2eV 13.605698  
#define au2a 0.52917
#define Ry2K 157886.6
#define P2GPa 14711.12
#define hc   (1)
#define re   (1) 
#define mec2 (1)
#endif 

/*#define amu2M 2.133086*/
/*#define Me   (0.00054858990945*amu2M)*/

#define t2fs 1.0

#define fs2t (1./t2fs)
#define M2amu (1.0/amu2M)
#define a2au (1./au2a)
#define eV2Ry (1./Ry2eV)
#define GPa2P (1./P2GPa)
#define K2Ry (1./(Ry2K))
#define hbar  (2e15*hbar_over_Eh_MKS) 
#define Qa  -1.0 
#define Me  ((me_MKS/u_MKS)*amu2M) 
#define kcoulomb 2.0 

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
