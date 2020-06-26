#ifndef BIOCHARMMTOP_H
#define	BIOCHARMMTOP_H

typedef struct mass_parms_str
{
    char atmType[5];
    int  atmTypeID;
    char element[3];
    double mass;
} MASS_PARMS;

typedef struct tatom_parms_str
{
    int atmID;
    char atmName[5];  
    MASS_PARMS *atmTypePtr;
    double charge;
//For Transforming molecules
    int species; 
} TATOM_PARMS;


typedef struct tgroup_parms_str
{
    int grpID;
    int grpAtomSize;
    TATOM_PARMS ** grpAtoms; 
    
} TGROUP_PARMS;

//typedef struct tone_parms_str
//{
//    TATOM_PARMS* atomI;      
//} TONE_PARMS;

//typedef struct ttwo_parms_str
//{
//    char atomI[5];
//    char atomJ[5];
//    //TATOM_PARMS* atomI;  
//    //TATOM_PARMS* atomJ;    
//} TTWO_PARMS;
//
//typedef struct tfour_parms_str
//{
//    char atomI[5];
//    char atomJ[5];
//    char atomK[5];
//    char atomL[5];    
//    //TATOM_PARMS* atomI;  
//    //TATOM_PARMS* atomJ;     
//    //TATOM_PARMS* atomK;  
//    //TATOM_PARMS* atomL;     
//   
//} TFOUR_PARMS;

typedef struct one_atom_str
{
    char atomI[5];     
} ONE_ATOM;

typedef struct two_atom_str
{
    char atomI[5];  
    char atomJ[5];    
} TWO_ATOM;

typedef struct four_atom_str
{
    char atomI[5];  
    char atomJ[5];
    char atomK[5];  
    char atomL[5];  
   
} FOUR_ATOM;

// I J K L R(I(J/K)) T(I(JK/KJ)) PHI T(JKL) R(KL)
typedef struct tic_parms_str
{
    FOUR_ATOM* atmsPtr;
    double kconst1;
    double angle1;
    double torsion;
    double angle2;
    double kconst2;
        
} TIC_PARMS;

typedef struct resi_parms_str
{
    int  resID;
    int  resType;
    char resName[5];    
    double charge;
    
    int atomListSize;
    int groupListSize;
    int bondListSize;
    int imprListSize;
    int cmapListSize;
    int donorListSize;
    int accepListSize;
    int delAtmListSize;
    int icListSize;
    int species1Size;
    int species2Size; 

    TATOM_PARMS ** atomList;
    TGROUP_PARMS ** groupList;        
    TWO_ATOM ** bondList;
    FOUR_ATOM ** imprList;
    FOUR_ATOM ** cmapList;
    ONE_ATOM ** donorList;
    ONE_ATOM ** accepList;
    ONE_ATOM ** delAtmList;
    ONE_ATOM ** species1List;
    ONE_ATOM ** species2List;
    TIC_PARMS ** icList;
    
} RESI_PARMS;

  
typedef struct charmm_top_str{
    int massParmSize;
    int resiParmSize;
    int presParmSize;

    MASS_PARMS ** massParms;
    
    RESI_PARMS ** resiParms;
    
    RESI_PARMS ** presParms;
    
}CHARMM_TOP;


int parseCharmmTop(const char *fileName, CHARMM_TOP *charmmTops);

int freeCharmmTop(CHARMM_TOP *charmmTops);

RESI_PARMS* findResiParms(CHARMM_TOP* charmmTop, char* resName);

double getMassbyAtmType(CHARMM_TOP* charmmTop, char* atmType);

#endif	/* BIOCHARMMTOP_H */

