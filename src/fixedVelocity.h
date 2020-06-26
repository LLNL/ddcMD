#ifndef FIXEDVELOCITY_H
#define FIXEDVELOCITY_H
#include "gid.h" 
#include "group.h" 
#include "collection.h" 
void fixedVelocity_parse(GROUP*g, gid_type label, char *field, int i);
char *fixedVelocity_write(GROUP*g, int i);
void fixedVelocity_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in);
void fixedVelocity_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt);
void fixedVelocity_parse(GROUP*g, gid_type label, char *field, int i);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
