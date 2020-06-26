#ifndef PAIR_ITERATOR_H
#define PAIR_ITERATOR_H

struct PairIterator_st;

typedef void (*advanceFcn) (struct PairIterator_st* iter);

typedef struct PairIterator_st
{
   int    ii;
   int    jj;
   double rx;
   double ry;
   double rz;
   double r2;

   advanceFcn advance;
   void* pairFinderState;
} PAIR_ITERATOR;

void pairIter_advance(PAIR_ITERATOR* iter);

#endif // #ifndef PAIR_ITERATOR_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
