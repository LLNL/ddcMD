#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <stdio.h> 
enum INTEGRATOR_CLASS { NGLF, NGLFNEW, NGLFNK,  NGLFRATTLE, NGLFCONSTRAINT, PNGLF, NGLFTEST, NVEGLF, NVEGLF_SIMPLE, NVTGLF, NPTGLF, STATIC, NEXTFILE, HYCOPINTEGRATOR };
typedef struct integrator_st
{
	char *name;
	char *objclass;
	char *value;
	char *type;
	void  *parent; 
	enum INTEGRATOR_CLASS itype;
	int uses_gpu;
	void (*eval_integrator) (void *, void *, void *parm);
	void (*writedynamic) (struct integrator_st *integrator, FILE *file);
	void *parms;
} INTEGRATOR;
INTEGRATOR *integrator_init(void *parent, char *name);
void integrator_writedynamic(INTEGRATOR*integrator, FILE*file);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
