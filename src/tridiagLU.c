#include <math.h>
#include <stdio.h>
void tridiagLU(int n, double *a, double *b, double *c, double *l, double *d, double *u) 
{
   d[0] = b[0]; 
   u[0] = c[0]; 
   l[0] = 0; 
   for (int i = 1;i<n-1;i++) 
   {
      l[i] = a[i]/d[i-1]; 
      d[i] = b[i] - l[i]*u[i-1]; 
      u[i] = c[i]; 
   }
   l[n-1] = a[n-1]/d[n-2]; 
   d[n-1] = b[n-1]-l[n-1]*u[n-2]; 
   u[n-1]=0; 
}
void tridiagLUSolve(int n, double *l, double *d, double *u,double *q, double *x)
{
   double y[n]; 
   y[0] =  q[0];
   for (int i=1;i<n;i++) 
   {
      y[i] = q[i] - l[i]*y[i-1]; 
   }
   x[n-1] = y[n-1]/d[n-1]; 
   for (int i=n-2;i>=0;i--)
   {
      x[i] = (y[i] - u[i]*x[i+1])/d[i]; 
   }
}
void tridiagSolveP(int n, double *a, double *b, double *c,double *r, double *x)
{
   double l[n],d[n],u[n]; 
   double x1[n],x2[n],r2[n]; 
   r2[0] = -a[0]; 
   for (int i=1;i<n-2;i++)  r2[i] =0; 
   r2[n-2] = -c[n-2]; 
   tridiagLU(n-1,a,b,c,l,d,u); 
   tridiagLUSolve(n-1, l, d, u, r, x1);
   tridiagLUSolve(n-1, l, d, u, r2, x2);

   double y= (r[n-1]-c[n-1]*x1[0]-a[n-1]*x1[n-2])/(b[n-1]+c[n-1]*x2[0]+a[n-1]*x2[n-2]); 

   for (int i=0;i<n-1;i++) x[i] =  x1[i]+y*x2[i]; 
   x[n-1] = y;
}

double  tridiagSolveCheck(int n, double *a, double *b, double *c, double *r, double *x) 
{
   double rS[n]; 
   
   int i = 0; 
   rS[i] = b[i] * x[i] + c[i] *x[i+1] ;
   for (i=1;i<n-1;i++) rS[i] = a[i]*x[i-1] + b[i] * x[i] + c[i] *x[i+1] ;
   i = n-1; 
   rS[i] = a[i]*x[i-1] + b[i] * x[i]  ;
   double sum = 0.0; 
   for (i=0;i<n;i++) 
   {
      double d= rS[i]-r[i]; 
      sum+=d*d; 
   }
   return sqrt(sum); 
}
double tridiagSolvePCheck(int n, double *a, double *b, double *c, double *r, double *x) 
{
   double rS[n]; 
   
   int i = 0; 
   rS[i] = a[i] *x[n-1] + b[i] * x[i] + c[i] *x[i+1] ;
   for (i=1;i<n-1;i++) rS[i] = a[i]*x[i-1] + b[i] * x[i] + c[i] *x[i+1] ;
   i = n-1; 
   rS[i] = a[i]*x[i-1] + b[i] * x[i]  + c[i]*x[0];
   double sum = 0.0; 
   for (i=0;i<n;i++) 
   {
      double d= rS[i]-r[i]; 
      sum+=d*d; 
   }
   return sqrt(sum); 
}
#ifdef SAT
#include <stdio.h>
#include <stdlib.h>
void matrixVector(int n, double **A, double x[], double r[])
{
   for (int i=0;i<n;i++)
   {
      double sum =0.0; 
      for (int j=0;j<n;j++)
      {
         sum += A[i][j]*x[j]; 
      }
      r[i] = sum; 
   }
}
void printMatix(int n, double **A)
{
   for (int i=0;i<n;i++)
   {
      double sum =0.0; 
      for (int j=0;j<n;j++)
      {
         printf(" %10.5f", A[i][j]); 
      }
      printf("\n"); 
   }
}
int main()
{
   int n = 100; 
   double a[n],b[n],c[n],q[n],r[n]; 
   double l[n],d[n],u[n]; 
   a[0] = -1; 
   b[0] = 2.1; 
   c[0] = -1; 
   q[0]  =  1.0;
   for (int i=1;i<n-1;i++)
   {
      a[i] = -1; 
      b[i] =  2.1; 
      c[i] = -1; 
      q[i]  = 1.0;  
   }
   a[n-1] = -1;  
   b[n-1] = 2.2; 
   c[n-1] =-1; 
   q[n-1]  = 1.0;        
   double x[n],x1[n],x2[n]; 

   tridiagLU(n-1,a,b,c,l,d,u); 

   r[0] = -a[0]; for (int i=1;i<n-2;i++)  r[i] =0; r[n-2] = -c[n-2]; 
   tridiagLUSolve(n-1, l, d, u, q, x1);
   tridiagLUSolve(n-1, l, d, u, r, x2);
   x[n-1] =  (q[n-1]-c[n-1]*x1[0]-a[n-1]*x1[n-2])/(b[n-1]+c[n-1]*x2[0]+a[n-1]*x2[n-2]); 
   for (int i=0;i<n-1;i++) x[i] = x1[i] + x[n-1]*x2[i]; 
   double res; 
   res=tridiagSolvePCheck(n, a, b, c, q, x) ;

   tridiagSolveP(n, a, b, c, q, x);
   res=tridiagSolvePCheck(n, a, b, c, q, x) ;


}
#endif
