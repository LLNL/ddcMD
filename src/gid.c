#include "gid.h"
int compareGid(const void *v1, const void *v2)
{
   if(*(gid_type*)v1 < *(gid_type*)v2) return -1;
   if(*(gid_type*)v1 > *(gid_type*)v2) return 1;
   return 0;
}

