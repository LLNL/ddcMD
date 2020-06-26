#ifndef FUNCTIONS_H
#define FUNCTIONS_H
enum FUNCTION_CLASS { LINEAR, STEPDOWN };
typedef struct function_st
{
	char *name;
	char *objclass;
	char *value;
	char *type;
	void *parent; 
	enum FUNCTION_CLASS itype;
	double (*f) ();
	void *parms;
} FUNCTION;


typedef struct linear_parm_st
{
	double m;
	double b;
} LINEAR_PARMS;

typedef struct stepdown_parm_st
{
	double x1, x2;
} STEPDOWN_PARMS;
#endif
