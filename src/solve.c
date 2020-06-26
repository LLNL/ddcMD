#include <math.h>
static void pivot(int n, double a[n][n], int indx[n]);
void solve(int n, double a[n][n],double *b,double *x)
{
  int indx[n];
  int i,j;

  pivot (n,a,indx);

  for(i = 0; i < n-1; ++i)
  {
    for(j = i+1; j < n; ++j)
    {
      b[indx[j]] = b[indx[j]]-a[indx[j]][i]*b[indx[i]];
    }
  }

  x[n-1] = b[indx[n-1]]/a[indx[n-1]][n-1];
  for (i = n-2; i>=0; i--)
  {
    x[i] = b[indx[i]];
    for (j = i+1; j < n; ++j)
    {
      x[i] = x[i]-a[indx[i]][j]*x[j];
    }
    x[i] = x[i]/a[indx[i]][i];
  }
}

void pivot (int n, double a[n][n], int indx[n])
{
  int i, j, k,itmp;
  double c1, pi, pi1, pj;
  double c[n];

/* Initialize the index */

  for (i = 0; i < n; ++i)
  {
    indx[i] = i;
  }

/* Find the rescaling factors, one from each row */
 
  for (i = 0; i < n; ++i)
  {
    c1 = 0;
    for (j = 0; j < n; ++j)
    {
      if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
    }
    c[i] = c1;
  }

/* Search the pivoting (largest) element from each column */ 

  k=0;
  for (j = 0; j < n-1; ++j)
  {
    pi1 = 0.0;
    for (i = j; i < n; ++i)
    {
      pi = fabs(a[indx[i]][j])/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

/* Interchange the rows via indx[] to record pivoting order */

    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < n; ++i)
    {
      pj = a[indx[i]][j]/a[indx[j]][j];

/* Record pivoting ratios below the diagonal */

      a[indx[i]][j] = pj;

/* Modify other elements accordingly */

      for (k = j+1; k < n; ++k)
      {
        a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
      }
    }
  }
}
/*
#include <stdio.h>
int main()
{
   double m[2][2];
   double b[]={1.0,1.0}; 
   double x[2]; 
   m[0][0]=1.0; 
   m[0][1]=0.5; 
   m[1][0]=0.25; 
   m[1][1]=1.5; 
   solve(2,m,b,x); 
   printf("%e %e\n",x[0],x[1]); 
}
int mainx()
{
	double  rc , aE,aD,aP,cP,c0[64*64],c1[64*64];
	int n=24; 
	rc = 1.957;
	aE = 1/.7268;
	aP = 1/0.06144;
	aD = sqrt(M_PI)/2.734e-4;
	cP =0; 
	double dvdrOverR; 
	
	plasma_DDBMakeSeries( n, rc, aE, aD, aP, cP, c0);
	double u0 = 1e-16;
 	double u1 =pow(0.5,nInterval-1); 
	for (double u=0.0;u<1.0;u+=1e-3)
	{
    	double v= plasma_DDBpair(c0,u*rc*rc,n, &dvdrOverR);
	printf("%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",sqrt(u)*rc,sqrt(u),u,v,fn(u,0),v.f-fn(u,0));
	}
	printf("end_of_data\n"); 
	aD = sqrt(M_PI)/7.379e-2;
	plasma_DDBMakeSeries( n, rc, aE, aD, aP, cP, c1);
	for (double u=0.0;u<1.0;u+=1e-3)
	{
    double v= plasma_DDBpair(c1,u*rc*rc,  n, &dvdrOverR);
	printf("%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",sqrt(u)*rc,sqrt(u),u,v,fn(u,0),v.f-fn(u,0));
	}
}
*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
