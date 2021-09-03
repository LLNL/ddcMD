#include "mpiTypes.h"

#include "domain.h"
#include "ddc.h" // PARTICLE

MPI_Datatype domainx_MPIType(void)
{
   static int intialized = 0;
   static MPI_Datatype domainxType;
   if (intialized == 0)
   {
      intialized = 1;
      MPI_Type_contiguous(sizeof(DOMAINX), MPI_BYTE, &domainxType);
      MPI_Type_commit(&domainxType);
   }
   return domainxType;
}
/*
MPI_Datatype domainx_MPIType(void)
{
   static int intialized = 0;
   static MPI_Datatype domainxType;
   if (intialized == 0)
   {
      intialized = 1;
      DOMAINX domain; 
      int n = 4; 
      int blkcnt[n]; 
      MPI_Aint disp[n]; 
      blkcnt[0]=1; 
      blkcnt[1]=1; 
      blkcnt[2]=1; 
      blkcnt[3]=1; 
      MPI_Datatype types[n]; 
      types[0] = MPI_DOUBLE; 
      types[1] = threeVector_MPIType();
      types[2] = threeVector_MPIType();
      types[3] = MPI_INT; 
      MPI_Get_address(&domain.radius, &disp[0]);
      MPI_Get_address(&domain.center, &disp[1]);
      MPI_Get_address(&domain.extent, &disp[2]);
      MPI_Get_address(&domain.shape , &disp[3]);
      for (int i = n-1; i >= 0; i--) disp[i] -= disp[0];
      MPI_Type_create_struct(n, blkcnt, disp, types, &domainxType);
      MPI_Type_commit(&domainxType);
   }
   return domainxType;
}
*/
MPI_Datatype threeVector_MPIType(void)
{
   static int intialized = 0;
   static MPI_Datatype threeVectorType;
   if (intialized == 0)
   {
      intialized = 1;
      MPI_Type_contiguous(3, MPI_DOUBLE, &threeVectorType);
      MPI_Type_commit(&threeVectorType);
   }
   return threeVectorType;
}
MPI_Datatype threeInt_MPIType(void)
{
   static int intialized = 0;
   static MPI_Datatype threeIntType;
   if (intialized == 0)
   {
      intialized = 1;
      MPI_Type_contiguous(3, MPI_INT, &threeIntType);
      MPI_Type_commit(&threeIntType);
   }
   return threeIntType;
}
MPI_Datatype particle_MPIType(void)
{
   static int intialized = 0;
   static MPI_Datatype particleType;
   if (intialized == 0)
   {
      intialized = 1;
      PARTICLE particle;
      int n = 3;
      int blkcnt[n];
      MPI_Aint disp[n];
      MPI_Datatype types[n];
      blkcnt[0] = 2;
      blkcnt[1] = 1;
      blkcnt[2] = 1;
      MPI_Get_address(&particle.domain_id, &disp[0]);
      MPI_Get_address(&particle.global_index, &disp[1]);
      MPI_Get_address(&particle.r.x, &disp[2]);
      types[0] = MPI_INT;
      types[1] = MPI_LONG_LONG;
      types[2] = threeVector_MPIType();
      for (int i = n-1; i >= 0; i--) disp[i] -= disp[0];
      MPI_Type_create_struct(n, blkcnt, disp, types, &particleType);
      MPI_Type_commit(&particleType);
   }
   return particleType;
}

   


/* Local Variables: */
/* tab-width: 3 */
/* End: */
