#ifndef COLLECTION_H
#define COLLECTION_H
#include "species.h"
#include "group.h"
#include "state.h"
#include "pio.h"
#include "gpunlist.h"

struct pfile_st; // forward declaration avoids including pio.h

enum COLLECTION_CLASS {
    DEFAULT
};

enum COLLECTION_WRITEMODE {
    NOWHEADER, WHEADER, WHEADER8FOLD, WREPLICATE
};

typedef struct collection_st {
    char *name; /* species name */
    char *objclass;
    char *value;
    char *type; /* model */
    void *parent;
    enum COLLECTION_CLASS itype; /* integer label for type */

    int size, max_imember;
    gid_type gsize;
    void *first, *last;    
    
    STATE *state;
    int dynamicInfoFromFileHeader;
    STATE *gpustate, *gpustate_h;
    GPUNLIST * gnlist;

} COLLECTION;

COLLECTION *collection_init(void *parent, char *name);

void collection_read(COLLECTION *collection, PFILE *pfile);
// private:

void collection_readASCII(STATE *state, int size, struct pfile_st *pfile);
void collection_readBINARY(STATE *state, int lengthRecord, int nRecords, struct pfile_st *pfile);
void collection_readFIXRECORDASCII(COLLECTION* c, struct pfile_st* pfile);
void collection_readVARRECORDASCII(COLLECTION* c, struct pfile_st* pfile);
void collection_readFIXRECORDBINARY(COLLECTION* c, struct pfile_st* pfile);


#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
