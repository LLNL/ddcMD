#ifndef CONSTRAINT_H
#define CONSTRAINT_H
typedef struct consPair_str
{
    int atomI;
    int atomJ;
    double d2;    
}CONSPAIR;

typedef struct constraintNew_str
{
    int numAtom;
    int numPair;
    int* atomList;
    CONSPAIR** consPair;
}CONSTRAINTNEW;
#endif
