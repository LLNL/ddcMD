#include "ddc.h"
static PARTICLE *p;

void sortsetparticle(PARTICLE*pvalue)
{
	p = pvalue;
}

int sortbyglobalindex(int *i1, int *i2)
{
	if (p[*i1].global_index > (p[*i2]).global_index) return 1;
	if (p[*i1].global_index < (p[*i2]).global_index) return -1;
	return 0;
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
