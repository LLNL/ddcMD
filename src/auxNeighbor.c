#include <stdlib.h>
#include "auxNeighbor.h"
#include "expandbuffer.h"
#include "error.h"
#include "object.h"
#include "mpiUtils.h"
static AUXNEIGHBOR *neighbor =NULL; 
AUXNEIGHBOR *auxNeighbor_provide(double rProvide)
{
	if (neighbor->rProvide > 0.0) 	return NULL; 
	if (rProvide< neighbor->rRequest) return NULL; 
	neighbor->rProvide = rProvide; 
	return neighbor; 
	
}
AUXNEIGHBOR *auxNeighbor_request(double rRequest)
{
	if (rRequest > neighbor->rRequest) neighbor->rRequest = rRequest;
	neighbor->numberPairs=0; 
	return neighbor; 
}
void auxNeighbor_begin(int n) 
{
	neighbor->numberPairs=0; 
	neighbor->rProvide=0.0; 
	neighbor->numberParticles=n; 
	neighbor->index=(AUXNEIGHBORINDEX *)ExpandBuffers((void *)neighbor->index, sizeof(AUXNEIGHBORINDEX), n, 1024, LOCATION("auxNeighbor_begin"), "neighbor->index");
	for (int i=0;i<n;i++) 
	{
		neighbor->index[i].startPairs=-1; 
		neighbor->index[i].nPairs=0; 
	}
}

AUXNEIGHBOR *auxNeighbor_expandList(int deltaNumber)
{
 	neighbor->list = (AUXNEIGHBORLIST *)ExpandBuffers((void *)neighbor->list, sizeof(AUXNEIGHBORLIST), neighbor->numberPairs+deltaNumber, 1024, LOCATION("auxNeighbor_expandList"), "neighbor->list");
	return neighbor; 
}
AUXNEIGHBOR *auxNeighbor_init(void *parent, char *name)
{
	object_compilestring("auxNeighbor AUXNEIGHBOR  {} " ) ;
	neighbor = (AUXNEIGHBOR*) object_initialize(name, "AUXNEIGHBOR", sizeof(AUXNEIGHBOR));
	neighbor->parent = NULL;  
	neighbor->rProvide = 0.0; 
	neighbor->rRequest = 0.0; 
	neighbor->numberPairs=0; 
	neighbor->list = NULL; 
	neighbor->index = NULL; 
	object_get((OBJECT *)neighbor, "rRequest", &neighbor->rRequest, WITH_UNITS, 1, "0","l",NULL);
	return neighbor; 
}
int auxNeighbor_sortByFirstIndex(AUXNEIGHBORLIST *a,AUXNEIGHBORLIST *b)
{
	if ( a->i < b->i ) return -1 ;
	if ( a->i > b->i ) return  1 ;
	if ( a->r2 < b->r2 ) return -1 ;
	if ( a->r2 > b->r2 ) return 1 ;
	return 0; 
}
AUXNEIGHBORLIST *auxNeighbor_list(void)
{
	qsort(neighbor->list,neighbor->numberPairs,sizeof(AUXNEIGHBORLIST),(int(*)(const void*,const void*))auxNeighbor_sortByFirstIndex);
	int ip =0; 
	int nPairs=0; 
	if (neighbor->numberPairs >  0 ) 
	{
		unsigned i = neighbor->list[ip].i; 
		neighbor->index[i].startPairs = ip; 
		nPairs++;
		for (ip=1;ip<neighbor->numberPairs;ip++) 
		{
			if ( neighbor->list[ip].i  != i )
			{
				neighbor->index[i].nPairs = nPairs; 
				nPairs=0; 
				i = neighbor->list[ip].i; 
				neighbor->index[i].startPairs = ip; 
			}
			nPairs++;
		}
		neighbor->index[i].nPairs = nPairs; 
	}
	return neighbor->list;
}
int auxNeighbor_numberPairs(void)
{
	return neighbor->numberPairs;
}
void auxNeighborFcn(AUXNEIGHBOR *neighbor,double r2, int ii,int jj) 
{
	if  ((neighbor) != NULL && (r2) < (neighbor->rRequest*neighbor->rRequest)) 
	{                                                
    	(neighbor)->list[neighbor->numberPairs].i =(ii);   
    	(neighbor)->list[neighbor->numberPairs].j =(jj); 
    	(neighbor)->list[neighbor->numberPairs].r2 =(r2); 
    	(neighbor)->list[neighbor->numberPairs+1].i =(jj);   
    	(neighbor)->list[neighbor->numberPairs+1].j =(ii); 
    	(neighbor)->list[neighbor->numberPairs+1].r2 =(r2); 
    	neighbor->numberPairs+=2;                       
	}
}

AUXNEIGHBORINDEX *auxNeighbor_index(void)
{
	return neighbor->index;
}
int auxNeighbor_numberParticles(void)
{
	return neighbor->numberParticles;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
