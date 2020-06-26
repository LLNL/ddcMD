#ifndef BIOMMFF_H
#define BIOMMFF_H

typedef struct massparms_st
{
    char *name;
    char *objclass;
    char *value;
    char *type;		/* model */
    void *parent;
    char *atomType;
    int atomTypeID;
    double mass;
    
} MASSPARMS;

MASSPARMS *massparms_init(void *parent, char *name);

typedef struct atomparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   int atomID;  
   char *atomName;
   char *atomType;
   int atomTypeID;
   double charge;
   
} ATOMPARMS;

ATOMPARMS *atomparms_init(void *parent, char *name);

typedef struct groupparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   int groupID; 
   
   int nAtoms;
   ATOMPARMS** atomList;
   
} GROUPPARMS;

GROUPPARMS *groupparms_init(void *parent, char *name);


typedef struct bondparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
 
   int atomI;
   int atomJ;
   int func;
   char * atomTypeI;
   char * atomTypeJ;
   double kb;
   double b0;
   
} BONDPARMS;

typedef struct consparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
 
   int atomI;
   int atomJ;
   char *atomTypeI;
   char *atomTypeJ;   
   int valid;   
   int func;
   double r0;
   
} CONSPARMS;


typedef struct conslistparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   
   int consSubListSize;
   CONSPARMS** consSubList;
   
}CONSLISTPARMS;

typedef struct excludeparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
 
   int atomI;
   int atomJ;
   char *atomTypeI;
   char *atomTypeJ;   
   int valid;
   
} EXCLUDEPARMS;

typedef struct angleparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
 
   int atomI;
   int atomJ;
   int atomK;
   int func;
   double ktheta;
   double theta0;
   
} ANGLEPARMS;


typedef struct torsparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
 
   int atomI;
   int atomJ;
   int atomK;
   int atomL;
   int func;
   int n;
   double kchi;
   double delta;
   
} TORSPARMS;

typedef struct resiparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   int resID;
   int resType;
   int centerAtom;
   char* resName;    
   double charge;
   int nGroups;
   int nBonds;
   int nCons;
   int nExclude;
   int nAngles;
   int nTors;
   GROUPPARMS** groupList;
   BONDPARMS** bondList;
   CONSLISTPARMS** consList;
   EXCLUDEPARMS** exclusionList;
   ANGLEPARMS** angleList;
   TORSPARMS** torsList;
   
} RESIPARMS;

RESIPARMS *resiparms_init(void *parent, char *name);

typedef struct ljparms_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent;  
   char *atomTypeI;
   char *atomTypeJ;
   int  indexI; 
   int  indexJ;
   double sigma;
   double eps;
   
} LJPARMS;

LJPARMS *ljparms_init(void *parent, char *name);

typedef struct mmff_st
{
   char *name;		
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   int nResiParms;
   int nLJParms;
   int nAtomType;
   
   RESIPARMS** resiParms;
   MASSPARMS** atomTypes;
   LJPARMS** ljParms;
   //void *parms;		/* pointer to  parameters for potential function */
} MMFF;

MMFF *mmff_init(void *parent,char *name);


#endif /* BIOMMFF_H */

