#include "bioMMFF.h"
#include "object.h"
#include "ddcMalloc.h"

MASSPARMS *massparms_init(void *parent, char *name)
{
    MASSPARMS *massparms;
    massparms = (MASSPARMS *) object_initialize(name, "MASSPARMS", sizeof (MASSPARMS));
    object_get((OBJECT *) massparms, "atomType", &(massparms->atomType), STRING, 1, "NoType");
    object_get((OBJECT *) massparms, "atomTypeID", &(massparms->atomTypeID), INT, 1, "0");
    object_get((OBJECT *) massparms, "mass", &(massparms->mass), WITH_UNITS, 1, "1.0", "m", NULL);
    return massparms;

}

ATOMPARMS *atomparms_init(void *parent, char *name)
{
    ATOMPARMS *atomparms;
    atomparms = (ATOMPARMS *) object_initialize(name, "ATOMPARMS", sizeof (ATOMPARMS));

    object_get((OBJECT *) atomparms, "atomID", &(atomparms->atomID), INT, 1, "0");
    object_get((OBJECT *) atomparms, "atomName", &(atomparms->atomName), STRING, 1, "NoName");
    object_get((OBJECT *) atomparms, "atomType", &(atomparms->atomType), STRING, 1, "NoType");
    object_get((OBJECT *) atomparms, "atomTypeID", &(atomparms->atomTypeID), INT, 1, "0");
    object_get((OBJECT *) atomparms, "charge", &(atomparms->charge), WITH_UNITS, 1, "0.0", "i*t", NULL);
    return atomparms;

}

GROUPPARMS *groupparms_init(void *parent, char *name)
{
    GROUPPARMS *groupparms;
    groupparms = (GROUPPARMS *) object_initialize(name, "GROUPPARMS", sizeof (GROUPPARMS));
    object_get((OBJECT *) groupparms, "groupID", &(groupparms->groupID), INT, 1, "0");

    char** atomNames;
    groupparms->nAtoms = object_getv((OBJECT *) groupparms, "atomList", (void *) &atomNames, STRING, ABORT_IF_NOT_FOUND);
    groupparms->atomList = (ATOMPARMS**) ddcMalloc(groupparms->nAtoms * sizeof (ATOMPARMS*));
    for (int i = 0; i < groupparms->nAtoms; i++)
    {
        //printf("%s\n", atomNames[i]);
        groupparms->atomList[i] = atomparms_init(groupparms, atomNames[i]);
    }

    return groupparms;
}

BONDPARMS *bondparms_init(void *parent, char *name)
{
    BONDPARMS *bondparms;
    bondparms = (BONDPARMS *) object_initialize(name, "BONDPARMS", sizeof (BONDPARMS));

    object_get((OBJECT *) bondparms, "atomI", &(bondparms->atomI), INT, 1, "0");
    object_get((OBJECT *) bondparms, "atomJ", &(bondparms->atomJ), INT, 1, "0");
    object_get((OBJECT *) bondparms, "func", &(bondparms->func), INT, 1, "1");
    object_get((OBJECT *) bondparms, "atomTypeI", &(bondparms->atomTypeI), STRING, 1, "NoType");
    object_get((OBJECT *) bondparms, "atomTypeJ", &(bondparms->atomTypeJ), STRING, 1, "NoType");
    object_get((OBJECT *) bondparms, "kb", &(bondparms->kb), WITH_UNITS, 1, "0.0", "kJ*mol^-1*nm^-2", NULL);
    object_get((OBJECT *) bondparms, "b0", &(bondparms->b0), WITH_UNITS, 1, "0.0", "nm", NULL);
    return bondparms;

}

CONSPARMS *consparms_init(void *parent, char *name)
{
    CONSPARMS *consparms;
    consparms = (CONSPARMS *) object_initialize(name, "CONSPARMS", sizeof (CONSPARMS));

    object_get((OBJECT *) consparms, "atomI", &(consparms->atomI), INT, 1, "0");
    object_get((OBJECT *) consparms, "atomJ", &(consparms->atomJ), INT, 1, "0");
    object_get((OBJECT *) consparms, "atomTypeI", &(consparms->atomTypeI), STRING, 1, "NoType");
    object_get((OBJECT *) consparms, "atomTypeJ", &(consparms->atomTypeJ), STRING, 1, "NoType");
    object_get((OBJECT *) consparms, "func", &(consparms->func), INT, 1, "1");
    object_get((OBJECT *) consparms, "r0", &(consparms->r0), WITH_UNITS, 1, "0.0", "nm", NULL);
    if (consparms->func == 1)
    {
        consparms->valid = 1;
    }
    else
    {
        consparms->valid = 0;
    }
    return consparms;

}

CONSLISTPARMS *conslistparms_init(void *parent, char *name)
{
    CONSLISTPARMS *conslistparms;
    conslistparms = (CONSLISTPARMS *) object_initialize(name, "CONSLISTPARMS", sizeof (CONSLISTPARMS));

    char** consNames;
    conslistparms->consSubListSize = object_getv((OBJECT *) conslistparms, "constraintSubList", (void *) &consNames, STRING, IGNORE_IF_NOT_FOUND); // bond is not required.
    conslistparms->consSubList = (CONSPARMS**) ddcMalloc(conslistparms->consSubListSize * sizeof (CONSPARMS*));
    for (int i = 0; i < conslistparms->consSubListSize; i++)
    {
        //printf("%s\n", bondNames[i]);
        conslistparms->consSubList[i] = consparms_init(conslistparms, consNames[i]);
    }

    return conslistparms;

}

EXCLUDEPARMS *excludeparms_init(void *parent, char *name)
{
    EXCLUDEPARMS *excludeparms;
    excludeparms = (EXCLUDEPARMS *) object_initialize(name, "EXCLUDEPARMS", sizeof (EXCLUDEPARMS));

    object_get((OBJECT *) excludeparms, "atomI", &(excludeparms->atomI), INT, 1, "0");
    object_get((OBJECT *) excludeparms, "atomJ", &(excludeparms->atomJ), INT, 1, "0");
    object_get((OBJECT *) excludeparms, "atomTypeI", &(excludeparms->atomTypeI), STRING, 1, "NoType");
    object_get((OBJECT *) excludeparms, "atomTypeJ", &(excludeparms->atomTypeJ), STRING, 1, "NoType");
    // Some exclusion pairs are duplicate in the bpairs. 
    // Assume we count all the pair here and turn off in the validateExclusions(MMFF *mmff) function if duplicated
    excludeparms->valid = 1;
    return excludeparms;

}

ANGLEPARMS *angleparms_init(void *parent, char *name)
{
    ANGLEPARMS *angleparms;
    angleparms = (ANGLEPARMS *) object_initialize(name, "ANGLEPARMS", sizeof (ANGLEPARMS));

    object_get((OBJECT *) angleparms, "atomI", &(angleparms->atomI), INT, 1, "0");
    object_get((OBJECT *) angleparms, "atomJ", &(angleparms->atomJ), INT, 1, "0");
    object_get((OBJECT *) angleparms, "atomK", &(angleparms->atomK), INT, 1, "0");
    object_get((OBJECT *) angleparms, "ktheta", &(angleparms->ktheta), WITH_UNITS, 1, "0.0", "kJ*mol^-1", NULL);
    object_get((OBJECT *) angleparms, "theta0", &(angleparms->theta0), DOUBLE, 1, "0");
    object_get((OBJECT *) angleparms, "func", &(angleparms->func), INT, 1, "1");
    return angleparms;

}

TORSPARMS *torsparms_init(void *parent, char *name)
{
    TORSPARMS *torsparms;
    torsparms = (TORSPARMS *) object_initialize(name, "TORSPARMS", sizeof (TORSPARMS));

    object_get((OBJECT *) torsparms, "atomI", &(torsparms->atomI), INT, 1, "0");
    object_get((OBJECT *) torsparms, "atomJ", &(torsparms->atomJ), INT, 1, "0");
    object_get((OBJECT *) torsparms, "atomK", &(torsparms->atomK), INT, 1, "0");
    object_get((OBJECT *) torsparms, "atomL", &(torsparms->atomL), INT, 1, "0");
    object_get((OBJECT *) torsparms, "func", &(torsparms->func), INT, 1, "1");
    object_get((OBJECT *) torsparms, "n", &(torsparms->n), INT, 1, "1");
    object_get((OBJECT *) torsparms, "kchi", &(torsparms->kchi), WITH_UNITS, 1, "0.0", "kJ*mol^-1", NULL);
    object_get((OBJECT *) torsparms, "delta", &(torsparms->delta), DOUBLE, 1, "0");
    return torsparms;

}

RESIPARMS *resiparms_init(void *parent, char *name)
{
    RESIPARMS *resiparms;
    resiparms = (RESIPARMS *) object_initialize(name, "RESIPARMS", sizeof (RESIPARMS));
    object_get((OBJECT *) resiparms, "resID", &(resiparms->resID), INT, 1, "0");
    object_get((OBJECT *) resiparms, "resType", &(resiparms->resType), INT, 1, "0");
    object_get((OBJECT *) resiparms, "centerAtom", &(resiparms->centerAtom), INT, 1, "0");
    object_get((OBJECT *) resiparms, "resName", &(resiparms->resName), STRING, 1, "NoName");
    object_get((OBJECT *) resiparms, "charge", &(resiparms->charge), WITH_UNITS, 1, "0.0", "i*t", NULL);

    //printf("%d %d %s %f\n", resiparms->resID, resiparms->resType, resiparms->resName, resiparms->charge);

    char** groupNames;
    resiparms->nGroups = object_getv((OBJECT *) resiparms, "groupList", (void *) &groupNames, STRING, ABORT_IF_NOT_FOUND);
    resiparms->groupList = (GROUPPARMS**) ddcMalloc(resiparms->nGroups * sizeof (GROUPPARMS*));
    for (int i = 0; i < resiparms->nGroups; i++)
    {
        //printf("%s\n", groupNames[i]);
        resiparms->groupList[i] = groupparms_init(resiparms, groupNames[i]);
    }

    char** bondNames;
    resiparms->nBonds = object_getv((OBJECT *) resiparms, "bondList", (void *) &bondNames, STRING, IGNORE_IF_NOT_FOUND); // bond is not required.
    resiparms->bondList = (BONDPARMS**) ddcMalloc(resiparms->nBonds * sizeof (BONDPARMS*));
    for (int i = 0; i < resiparms->nBonds; i++)
    {
        //printf("%s\n", bondNames[i]);
        resiparms->bondList[i] = bondparms_init(resiparms, bondNames[i]);
    }

    char** conslistNames;
    resiparms->nCons = object_getv((OBJECT *) resiparms, "constraintList", (void *) &conslistNames, STRING, IGNORE_IF_NOT_FOUND); // constraint is not required.
    resiparms->consList = (CONSLISTPARMS**) ddcMalloc(resiparms->nCons * sizeof (CONSLISTPARMS*));
    for (int i = 0; i < resiparms->nCons; i++)
    {
        //printf("%s\n", consNames[i]);
        resiparms->consList[i] = conslistparms_init(resiparms, conslistNames[i]);
    }

    char** excludeNames;
    resiparms->nExclude = object_getv((OBJECT *) resiparms, "exclusionList", (void *) &excludeNames, STRING, IGNORE_IF_NOT_FOUND); // exclusion is not required.
    resiparms->exclusionList = (EXCLUDEPARMS**) ddcMalloc(resiparms->nExclude * sizeof (EXCLUDEPARMS*));
    for (int i = 0; i < resiparms->nExclude; i++)
    {
        //printf("%s\n", consNames[i]);
        resiparms->exclusionList[i] = excludeparms_init(resiparms, excludeNames[i]);
    }

    char** angleNames;
    resiparms->nAngles = object_getv((OBJECT *) resiparms, "angleList", (void *) &angleNames, STRING, IGNORE_IF_NOT_FOUND); // bond is not required.
    resiparms->angleList = (ANGLEPARMS**) ddcMalloc(resiparms->nAngles * sizeof (ANGLEPARMS*));
    for (int i = 0; i < resiparms->nAngles; i++)
    {
        //printf("%s\n", angleNames[i]);
        resiparms->angleList[i] = angleparms_init(resiparms, angleNames[i]);
    }

    char** torsNames;
    resiparms->nTors = object_getv((OBJECT *) resiparms, "dihedralList", (void *) &torsNames, STRING, IGNORE_IF_NOT_FOUND); // bond is not required.
    resiparms->torsList = (TORSPARMS**) ddcMalloc(resiparms->nTors * sizeof (TORSPARMS*));
    for (int i = 0; i < resiparms->nTors; i++)
    {
        //printf("%s\n", torsNames[i]);
        resiparms->torsList[i] = torsparms_init(resiparms, torsNames[i]);
    }

    return resiparms;
}

LJPARMS *ljparms_init(void *parent, char *name)
{
    LJPARMS *ljparms;
    ljparms = (LJPARMS *) object_initialize(name, "LJPARMS", sizeof (LJPARMS));

    object_get((OBJECT *) ljparms, "atomtypeI", &(ljparms->atomTypeI), STRING, 1, "NoType");
    object_get((OBJECT *) ljparms, "atomtypeJ", &(ljparms->atomTypeJ), STRING, 1, "NoType");
    object_get((OBJECT *) ljparms, "indexI", &(ljparms->indexI), INT, 1, "0");
    object_get((OBJECT *) ljparms, "indexJ", &(ljparms->indexJ), INT, 1, "0");
    object_get((OBJECT *) ljparms, "sigma", &(ljparms->sigma), WITH_UNITS, 1, "1.0", "nm", NULL);
    object_get((OBJECT *) ljparms, "eps", &(ljparms->eps), WITH_UNITS, 1, "0.0", "kJ*mol^-1", NULL);
    return ljparms;
}

MMFF *mmff_init(void *parent, char *name)
{
    MMFF *mmff;
    //char *type;
    char** resiParmNames;

    mmff = (MMFF *) object_initialize(name, "MMFF", sizeof (MMFF));
    mmff->parent = parent;

    mmff->nResiParms = object_getv((OBJECT *) mmff, "resiParms", (void *) &resiParmNames, STRING, ABORT_IF_NOT_FOUND);
    mmff->resiParms = (RESIPARMS**) ddcMalloc(mmff->nResiParms * sizeof (RESIPARMS*));
    for (int i = 0; i < mmff->nResiParms; i++)
    {
        //printf("%s\n", resiParmNames[i]);
        mmff->resiParms[i] = resiparms_init(mmff, resiParmNames[i]);
    }

    char** atomTypeNames;
    mmff->nAtomType = object_getv((OBJECT *) mmff, "atomTypeList", (void *) &atomTypeNames, STRING, ABORT_IF_NOT_FOUND);
    //mmff->atomTypes = (MASSPARMS**) ddcMalloc(mmff->nLJParms * sizeof (MASSPARMS*));
    mmff->atomTypes = (MASSPARMS**) ddcMalloc(mmff->nAtomType * sizeof (MASSPARMS*));
    for (int i = 0; i < mmff->nAtomType; i++)
    {
        // printf("%s\n", atomTypeNames[i]);
        mmff->atomTypes[i] = massparms_init(mmff, atomTypeNames[i]);
    }

    char** ljParmNames;
    mmff->nLJParms = object_getv((OBJECT *) mmff, "ljParms", (void *) &ljParmNames, STRING, ABORT_IF_NOT_FOUND);
    mmff->ljParms = (LJPARMS**) ddcMalloc(mmff->nLJParms * sizeof (LJPARMS*));
    for (int i = 0; i < mmff->nLJParms; i++)
    {
        //printf("%s\n", ljParmNames[i]);
        mmff->ljParms[i] = ljparms_init(mmff, ljParmNames[i]);
    }

    return mmff;
}
