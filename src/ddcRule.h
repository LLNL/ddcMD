#ifndef DDCRULE_H
#define DDCRULE_H

enum DDCRULE_CLASS
{
    DEFAULT_NONE, DDC_CHARMM, DDC_MARTINI, DDC_MOLECULE
};

typedef struct ddcrule_st
{
    char *name;
    char *objclass;
    char *value;
    char *type;
    void *parent;
    enum DDCRULE_CLASS itype;
    double firstRcut;
    void (*eval_ddcRule)(void* ddc);
    void *parms;

} DDCRULE;

DDCRULE *ddcRule_init(void *parent, char* name);

#endif /* DDCRULE_H */

