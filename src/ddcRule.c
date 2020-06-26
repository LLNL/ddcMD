#include <string.h>

#include "ddcRule.h"
#include "object.h"
#include "ddcMalloc.h"
#include "ddc.h"

#include "simulate.h"
#include "utilities.h"

#include "bioCharmmRule.h"
#include "bioMartiniRule.h"
#include "ddcRuleMolecule.h"

void ddcRuleDefault(DDC* ddc);

DDCRULE *ddcRule_init(void *parent, char* name)
{
    DDCRULE *ddcRule;
    ddcRule = (DDCRULE *) object_initialize(name, "DDCRULE", sizeof (DDCRULE));
    ddcRule->parent = parent;

    if (strcmp(name, "defaultRue") == 0)
    { //if it DDCRULE is missing in the object.data
        strcpy(ddcRule->type, "DEFAULT");
    }

    char* type;
    object_get((OBJECT *) ddcRule, "type", &type, STRING, 1, "DEFAULT");
    object_get((OBJECT *) ddcRule, "firstRcut", &(ddcRule->firstRcut), WITH_UNITS, 1, "40.0", "Angstrom", NULL);

    ddcRule->type = strdup(type);

    if (strcmp(type, "CHARMM") == 0)
    {
        ddcRule->itype = DDC_CHARMM;
        ddcRule->eval_ddcRule = (void (*)(void*))ddcRuleCharmm;
        ddcRuleCharmm_init(ddcRule);

    }
    else if (strcmp(type, "MARTINI") == 0)
    {
        ddcRule->itype = DDC_MARTINI;
        ddcRule->eval_ddcRule = (void (*)(void*))ddcRuleMartini;
        ddcRuleMartini_init(ddcRule);

    }
    else if (strcmp(type, "MOLECULE") == 0)
    {
        ddcRule->itype = DDC_MOLECULE;
        ddcRule->eval_ddcRule = (void (*)(void*))ddcRuleMolecule;
        ddcRuleMolecule_init(ddcRule);

    }
    else
    { //Default
        ddcRule->itype = DEFAULT_NONE;
        ddcRule->eval_ddcRule = (void (*)(void*))ddcRuleDefault;
        ddcRule->parms = NULL;
    }
    //ddcFree(type);
    return ddcRule;
}

// Default - do nothing

void ddcRuleDefault(DDC* ddc)
{

}

