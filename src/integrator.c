#include "integrator.h"
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include "three_algebra.h"
#include "object.h"
#include "nglf.h"
#include "nglfNew.h"
#include "nglfNK.h"
#include "hycopIntegrator.h"
#include "nglfrattle.h"
#include "nglfconstraint.h"
#include "nptglf.h"
#include "pnglf.h"
#include "ddcMalloc.h"
#include "HAVEGPU.h"

//void *nglfconstraint_parms(INTEGRATOR*integrator);

void *nptglf_parms(INTEGRATOR*integrator);
void *nveglf_parms(INTEGRATOR*integrator);
void *nveglf_simple_parms(INTEGRATOR*integrator);
void *static_parms(INTEGRATOR*integrator);
void *nextfile_parms(INTEGRATOR*integrator);
void nptglf(void *, void *simulate, void *p);
void nveglf(void *, void *simulate, void *p);
//void nglfconstraint(void *, void *simulate, void *p);
void nveglf_simple(void *, void *simulate, void *p);
void nextfile(void *, void *simulate, void *p);

void nptglf_writedynamic(INTEGRATOR *integrator,FILE*file);
static void writeNULL(INTEGRATOR *integrator,FILE*file)
{
}
INTEGRATOR *integrator_init(void *parent, char *name)
{
	char *type;
	static INTEGRATOR *integrator;

	integrator = (INTEGRATOR *) object_initialize(name, "INTEGRATOR", sizeof(INTEGRATOR));
	integrator->parent = parent; 
	object_get((OBJECT *) integrator, "type", &type, STRING, 1, NULL);

	integrator->name = ddcCalloc(strlen(name) + 1, sizeof(char));
	integrator->type = ddcCalloc(strlen(type) + 1, sizeof(char));
	integrator->writedynamic = writeNULL; 
	strcpy(integrator->name, name);
	strcpy(integrator->type, type);
   integrator->uses_gpu=0;
	if (strcmp(type, "NPTGLF") == 0)
	{
		integrator->itype = NPTGLF;
		integrator->parms = (void *)nptglf_parms(integrator);
		integrator->eval_integrator = (void (*)(void *, void*, void *parm)) nptglf;
		integrator->writedynamic = nptglf_writedynamic; 
	}
	if (strcmp(type, "NVTGLF") == 0 || strcmp(type, "NGLF") == 0 ) 
	{
		integrator->itype = NGLF;
		integrator->parms = (void *)nglf_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglf;
	}
	if (strcmp(type, "NGLFNEW") == 0 ) 
	{
		integrator->itype = NGLFNEW;
		integrator->parms = (void *)nglfNew_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfNew;
	}
	if (strcmp(type, "NGLFNK") == 0 ) 
	{
		integrator->itype = NGLFNK;
		integrator->parms = (void *)nglfNK_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfNK;
	}
	if (strcmp(type, "NGLFGPU") == 0 ) 
	{
		integrator->itype = NGLF;
		integrator->uses_gpu = 1;
		integrator->parms = (void *)nglf_parms(integrator);
		GPUCODE(integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfGPU;)
	}
	if (strcmp(type, "NGLFGPULANGEVIN") == 0 ) 
	{
		integrator->itype = NGLF;
		integrator->uses_gpu = 1;
		integrator->parms = (void *)nglf_parms(integrator);
		GPUCODE(integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfGPULangevin;)
	}
	if (strcmp(type, "HYCOP") == 0 ) 
	{
		integrator->itype = HYCOPINTEGRATOR;
		integrator->parms = (void *)hycopIntegrator_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))hycopIntegrator;
	}
	if (strcmp(type, "NGLFRATTLE") == 0) 
	{
		integrator->itype = NGLFRATTLE;
		integrator->parms = (void *)nglfrattle_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfrattle;
	} 
	if (strcmp(type, "NGLFCONSTRAINT") == 0) 
	{
		integrator->itype = NGLFCONSTRAINT;
		integrator->parms = (void *)nglfconstraint_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfconstraint;
	}        
	if (strcmp(type, "NGLFCONSTRAINTGPU") == 0) 
	{
		integrator->itype = NGLFCONSTRAINT;
      integrator->uses_gpu = 1;
		integrator->parms = (void *)nglfconstraint_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfconstraintGPU;
	} 
	if (strcmp(type, "NGLFCONSTRAINTGPULANGEVIN") == 0 ) 
	{
		integrator->itype = NGLFCONSTRAINT;
		integrator->uses_gpu = 1;
		integrator->parms = (void *)nglfconstraint_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfconstraintGPULangevin;
	}
	if (strcmp(type, "NGLFCONSTRAINTGPULANGEVINLCG64") == 0 ) 
	{
		integrator->itype = NGLFCONSTRAINT;
		integrator->uses_gpu = 1;
		integrator->parms = (void *)nglfconstraint_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfconstraintGPULangevinLCG64;
	}
	if (strcmp(type, "PNGLF") == 0 ) 
	{
		integrator->itype = PNGLF;
		integrator->parms = (void *)pnglf_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))pnglf;
	}
	if (strcmp(type, "NGLFERROR") == 0 ) 
	{
		integrator->itype = NGLF;
		integrator->parms = (void *)nglfError_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfError;
	}
	if (strcmp(type, "NGLFTEST") == 0 ) 
	{
		integrator->itype = NGLFTEST;
		integrator->parms = (void *)nglfTest_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nglfTest;
	}
	if (strcmp(type, "NVEGLF") == 0)
	{
		integrator->itype = NVEGLF;
		integrator->parms = (void *)nveglf_parms(integrator);
		integrator->eval_integrator = (void (*) (void *, void *, void *parm))nveglf;
	}
	if (strcmp(type, "NVEGLF_SIMPLE") == 0)
	{
		integrator->itype = NVEGLF_SIMPLE;
		integrator->parms = (void *)nveglf_simple_parms(integrator);
		integrator->eval_integrator = (void (*)(void *, void*, void *parm)) nveglf_simple;
	}
	if (strcmp(type, "NEXTFILE") == 0)
	{
		integrator->itype = NEXTFILE;
		integrator->eval_integrator = (void (*)(void *, void*, void *parm)) nextfile;
		integrator->parms = (void *)nextfile_parms(integrator);
	}
	return integrator;
}

void integrator_writedynamic(INTEGRATOR*integrator, FILE*file)
{
	switch (integrator->itype)
	{
	case NPTGLF:
		fprintf(file, "%s INTEGRATOR { zeta=%16.12e ; }\n", integrator->name, ((NPTGLF_PARMS *) integrator->parms)->zeta);
		break;
	default:
		break ;
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
