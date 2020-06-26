#ifndef AUXNEIGHBOR_H
#define AUXNEIGHBOR_H
typedef struct auxneighborIndex_str { int startPairs,nPairs;} AUXNEIGHBORINDEX ; 
typedef struct auxneighborList_str { unsigned int i,j; double r2;} AUXNEIGHBORLIST;
typedef struct auxneighbor_str
{ 
	char *name;
	char *objclass;
	char *value;
	char *type;		
	void  *parent; 
	double rRequest; 
	double rProvide; 
	int numberPairs;
	int numberParticles;
	AUXNEIGHBORLIST *list;
	AUXNEIGHBORINDEX *index; 
} AUXNEIGHBOR;
#define AddNeighbor(neighbor,r2,ii,jj)					 \
if  ((neighbor) != NULL && (r2) < (neighbor->rRequest*neighbor->rRequest)) \
{  										     	 \
	(neighbor)->list[neighbor->numberPairs].i =(ii);   \
	(neighbor)->list[neighbor->numberPairs].j =(jj); \
	(neighbor)->list[neighbor->numberPairs].r2 =(r2); \
	(neighbor)->list[neighbor->numberPairs+1].i =(jj);   \
	(neighbor)->list[neighbor->numberPairs+1].j =(ii); \
	(neighbor)->list[neighbor->numberPairs+1].r2 =(r2); \
	neighbor->numberPairs+=2;						\
}
AUXNEIGHBOR *auxNeighbor_init(void *parent, char *name);
void auxNeighbor_begin(int);
AUXNEIGHBOR *auxNeighbor_provide(double rProvide);
AUXNEIGHBOR *auxNeighbor_request(double rRequest);
AUXNEIGHBOR *auxNeighbor_expandList(int deltaNumber);
AUXNEIGHBORLIST *auxNeighbor_list(void);
int auxNeighbor_numberPairs(void);
int auxNeighbor_numberParticles(void);
AUXNEIGHBORINDEX *auxNeighbor_index(void);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
