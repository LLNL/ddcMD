#ifndef EAM_SC_H
#define EAM_SC_H

/** Sutton-Chen form
 *
 * \varphi(r_{ij}) = \epsilon(a/r_{ij})^n
 *
 * \rho(r_{ij}) = (a/r_{ij})^m
 *
 * F(\rho) = -c\epsilon\sqrt{\rho}
 *
 * Parameters are epsilon (eV), a (Angstrom), n, m, and c.  We also have
 * the ususal eam rmax parameter that is read and stored by the "base
 * class".
*/

#include "object.h"
#include "potential.h"
#include "eam.h"

void eam_sc_parms(POTENTIAL* object, EAM_PARMS* parms);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
