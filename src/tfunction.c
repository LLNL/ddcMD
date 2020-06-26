#include <stdio.h>
#include <math.h>
#include "tfunction.h"
#include "spline.h"
#include "object.h"
#include "mpiUtils.h"
#include "units.h"
#include "ddcMalloc.h"
#include "external.h"

//static void tf_objectNotFoundError(const char* name);

static TFUNCTION tfunc_read(const char* name);

#ifdef WITH_MPI
void tfunc_bcast(TFUNCTION* tab, int root, MPI_Comm comm);
#endif


TFUNCTION tfunc_init(const char* name)
{
   //OBJECT* obj = object_find(name, "TFUNCTION");
   //if (obj == NULL)  tf_objectNotFoundError(name);

   TFUNCTION table;
   if (getRank(0) == 0) 
      table = tfunc_read(name);
   #ifdef WITH_MPI
   tfunc_bcast(&table, 0, COMM_LOCAL);
   #endif
   return table;
}

#ifdef WITH_MPI
void tfunc_bcast(TFUNCTION* tab, int root, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   MPI_Bcast(&tab->type, 1, MPI_INT, root, comm);
   MPI_Bcast(&tab->rows, 1, MPI_INT, root, comm);
   MPI_Bcast(&tab->cols, 1, MPI_INT, root, comm);
   MPI_Bcast(&tab->x0, 1, MPI_DOUBLE, root, comm);
   MPI_Bcast(&tab->dx, 1, MPI_DOUBLE, root, comm);
  
   if (tab->type==1)
   {
      if (myRank != root)
      {
         tab->coef0 = ddcMalloc(4*(tab->rows-1)*sizeof(double));
      }
      MPI_Bcast(tab->coef0, 4*(tab->rows-1), MPI_DOUBLE, root, comm);
   }
   if (tab->type==2)
   {
      if (myRank != root)
      {
         tab->coef1 = ddcMalloc(4*(tab->rows-1)*sizeof(double));
         tab->coef2 = ddcMalloc(4*(tab->rows-1)*sizeof(double));
      }
      MPI_Bcast(tab->coef1, 4*(tab->rows-1), MPI_DOUBLE, root, comm);
      MPI_Bcast(tab->coef2, 4*(tab->rows-1), MPI_DOUBLE, root, comm);
   }
}
#endif

void tfunc_f(TFUNCTION* tab, double x, double* f1, double* f2)
{
   int i;
   double xi, X;
   X=((x-tab->x0)/tab->dx);
   i=(int)(X);
   xi=tab->x0+i*tab->dx;
   X=(x-xi);
   *f1=((*(tab->coef1+4*i)*X + *(tab->coef1+4*i+1))*X + *(tab->coef1+4*i+2))*X + *(tab->coef1+4*i+3);
   *f2=((*(tab->coef2+4*i)*X + *(tab->coef2+4*i+1))*X + *(tab->coef2+4*i+2))*X + *(tab->coef2+4*i+3);
}

void tfunc_df(TFUNCTION* tab, double x, double* df1, double* df2)
{
   int i;
   double xi, X;
   X=((x-tab->x0)/tab->dx);
   i=(int)(X);
   xi=tab->x0+i*tab->dx;
   X=(x-xi);
   *df1=2.*((*(tab->coef1+4*i)*3.*X + *(tab->coef1+4*i+1)*2.)*X + *(tab->coef1+4*i+2));
   *df2=2.*((*(tab->coef2+4*i)*3.*X + *(tab->coef2+4*i+1)*2.)*X + *(tab->coef2+4*i+2));
}

void tfunc_fdf(TFUNCTION* tab, double x, double* f, double* df)
{
   int i;
   double xi, X;
   X=((x-tab->x0)/tab->dx);
   i=(int)(X);
   xi=tab->x0+i*tab->dx;
   X=x-xi;
   *f=((*(tab->coef0+4*i)*X + *(tab->coef0+4*i+1))*X + *(tab->coef0+4*i+2))*X + *(tab->coef0+4*i+3);
   *df=(*(tab->coef0+4*i)*3.*X + *(tab->coef0+4*i+1)*2.)*X + *(tab->coef0+4*i+2);
}


//void tf_objectNotFoundError(const char* name)
//{
//   if (getRank(0) == 0)
//   {
//      printf("ERROR:  Can't find TFUNCTION object with name %s\n", name);
//      abortAll(-1);
//   }
//}


TFUNCTION tfunc_read(const char* name)
{
   TFUNCTION table= {0, 0, 0, 0.0, 0.0, NULL, NULL, NULL};
   
   unsigned type, rows, cols;
   double x0, dx;
   FILE *file = fopen(name, "r");
   if (file == NULL && getRank(0) == 0)
   {
      printf("FATAL ERROR:  Can't open table file %s\n", name);
      abortAll(-1);
   }
   fscanf(file, "%d",  &type);
   fscanf(file, "%d",  &rows);
   fscanf(file, "%d",  &cols);
   fscanf(file, "%lf", &x0);
   fscanf(file, "%lf", &dx);

   double scratch1[rows];
   double scratch2[rows];
   
   double length_convert = units_convert(1.0,"Angstrom",NULL); 
   double energy_convert = units_convert(1.0,"eV",NULL);
   if (type==1)//embedding energy
   {
      for (unsigned i=0;i<rows;++i)
      {
         fscanf(file, "%lf", &scratch1[i]);
         scratch1[i] *= energy_convert;
      }   
      table.coef0 = ddcMalloc(4*(rows-1) * sizeof(double));
      table.coef1 = NULL;
      table.coef2 = NULL;
      splcoef(dx*pow(length_convert, 2.0),scratch1,rows,table.coef0);
      table.type=type;
      table.rows=rows;
      table.cols=cols;
      table.x0=x0*pow(length_convert, 2.0);
      table.dx=dx*pow(length_convert, 2.0);
   }
   if (type==2)//pair potential and density
   {
      for (unsigned i=0;i<rows;++i)
      {
         fscanf(file, "%lf", &scratch1[i]);
         scratch1[i] *= energy_convert;
         fscanf(file, "%lf", &scratch2[i]);
         scratch2[i] *= pow(length_convert, 2.0);
      }   
      table.coef0 = NULL;
      table.coef1 = ddcMalloc(4*(rows-1) * sizeof(double));
      table.coef2 = ddcMalloc(4*(rows-1) * sizeof(double));
      splcoef(dx*pow(length_convert, 2.0),scratch1,rows,table.coef1);
      splcoef(dx*pow(length_convert, 2.0),scratch2,rows,table.coef2);
      table.type=type;
      table.rows=rows;
      table.cols=cols;
      table.x0=x0*pow(length_convert, 2.0);
      table.dx=dx*pow(length_convert, 2.0);
   }
   fclose(file);
   return table;
}
#ifdef SA
#include "codata.h"
int main(int argc, char *argv[])
{
   units_internal(a0_MKS,Rinfhc_MKS*1e-30/(a0_MKS*a0_MKS),1e-15,e_MKS/1e-15,Rinfhc_eV/kB_eV,1.0,1.0);
   units_external(1e-10,u_MKS,1e-15,e_MKS/1e-15,1.0,1.0,1.0);
   double cL2 = units_convert(1.0,NULL,"Angstrom^2"); 
   double cE = units_convert(1.0,NULL,"eV");
   char *name = "table.data";
   TFUNCTION tf = tfunc_init(name);
   double x0 = tf.x0;
   double x1 = x0 + tf.dx*(tf.rows-1);
   x1 = x0 + 11/cL2; 
   x0 += 10/cL2; 
   double dx = 0.01; 
   int n = (x1 - x0)/dx +1.000001; 

   double x = x0; 
   double s2  = (3.4*3.4)/(x1*cL2); 
   double s6 = s2*s2*s2; 
   double s12 = s6*s6; 
   double f0 = 4*1e-2*(s12-s6); 
   while (x <= x1*1.00001) 
   {
      double f1,f2;
      double df1,df2;
      tfunc_f(&tf, x, &f1, &f2);
      tfunc_df(&tf, x, &df1, &df2);
      double r2 = x *cL2; 
      double s2  = (3.4*3.4)/r2; 
      double s6 = s2*s2*s2; 
      double s12 = s6*s6; 
      f2 = 4e-2*(s12-s6)-f0; 
      df2 =-4e-2*(6*s12 - 3*s6)/r2;
      //df1 *= 0.5*cE/cL2; 
      df1 *=0.5;
      //printf("%f %e %e %e\n",r2,df1,df2,(df1-df2)/fabs(df2)); 
      printf("%f %e %e\n",x,f1,df1); 
      x += dx; 
   }
}
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
