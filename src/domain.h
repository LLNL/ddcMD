#ifndef DOMAIN_H
#define DOMAIN_H

#include "three_algebra.h"
#include <semaphore.h>

#ifdef __cplusplus
extern "C" {
#endif

   /** We have a circular coupling problem with most of the functions in
    * domain.c  Nearly all of these functions take a DDC* as an argument.
    * So, this file would need to include ddc.h.  However, the DDC
    * structure has a DOMAINX as a member so ddc.h needs to include this
    * file.  For now, to avoid the circular dependency we'll just avoid
    * declaring most of the domain related functions.  Eventually we may
    * refactor the domain functions to avoid the DDC* argument or come up
    * with some other solution.
    */
   enum DSHAPE {SPHERE, BRICK}; 

   typedef struct domain_st
   {
      double radius;
      THREE_VECTOR center;
      THREE_VECTOR extent;
      int shape; 
   }
   DOMAINX;

   typedef struct domainset_st
   {
      unsigned local_id;
      int      size;

      DOMAINX* domains;
      DOMAINX* local_domain; // If WITH_SHM is off, local_domain points into domains.

      /* Stuff not to be touched outside domain.c */
   }
   DOMAINSET;

   //unused
   //int domain_overlap(DOMAINX d1, DOMAINX d2, double rcut);
   int domain_possibleRemote(const DOMAINX d, THREE_VECTOR r, double rcut);
   double domain_separating_plane(THREE_VECTOR *c1, THREE_VECTOR *c2,THREE_VECTOR *u);
   // unused
   //void domain_print(DOMAINX domain);
   int domainset_nearestDomain(DOMAINSET* domainset, THREE_VECTOR r);
   THREE_VECTOR domainset_getLocalCenter(DOMAINSET* domainset);
   THREE_VECTOR* domainset_getPointerLocalCenter(DOMAINSET* domainset);
   THREE_VECTOR domainset_getRemoteCenter(DOMAINSET* domainset, int remote_id);
   void domainset_setLocalCenter(DOMAINSET* domainset, THREE_VECTOR new_center);
   double domainset_getLocalRadius(DOMAINSET* domainset);
   void domainset_setLocalRadius(DOMAINSET* domainset, double radius);
   void domainset_setLocalShape(DOMAINSET* domainset, int shape);
   THREE_VECTOR domainset_getLocalExtent(DOMAINSET* self);
   void domainset_setLocalExtent(DOMAINSET* self, THREE_VECTOR extent);
   void domainset_setLocalShape(DOMAINSET* self, int shape);

   void domainset_setRadius(DOMAINSET* domainset, int id, double radius);
   int domainset_overlap(DOMAINSET* domainset, int remote_id, double rcut);
   int domainset_overlapBox(DOMAINSET* domainset, int remote_id, double rcut);
   int domainset_overlapSphere(DOMAINSET* domainset, int remote_id, double rcut);
   void domainset_init(DOMAINSET* dset);
   int domain_side(DOMAINX* d1, DOMAINX* d2);
   void domainset_allGather(DOMAINSET* domainset);
   double domainset_getMaxRadius(DOMAINSET* domainset);
   double domainset_setRadiusToMax(DOMAINSET* domainset);
   double domainset_separatingPlane(DOMAINX* local_domain, DOMAINX *remote_domain, THREE_VECTOR *u);
   int domainset_particle(DOMAINSET* domainset, THREE_VECTOR v, int n, int *list);
   int domainset_particleAll(DOMAINSET* domainset, THREE_VECTOR v);
   void domainset_randomCenters(DOMAINSET* mydomains, int lx, int ly, int lz);
   void domain_getBoundingBox(THREE_VECTOR* rMin, THREE_VECTOR* rMax);

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
