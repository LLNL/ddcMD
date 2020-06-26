#include "spline.h"
#include <stdio.h>
#include <stdlib.h>

void splcoef(double dx, double y[], int rows, double coef[])
{
  int i;
  double cp[rows-2], dp[rows-2], M[rows];
  cp[0]=1./4.;
  dp[0]=(6./(dx*dx))*(y[0]-2.*y[1]+y[2])/4.;
  for (i=1;i<rows-2;i++)
  {
    cp[i]=1./(4.-cp[i-1]);
    dp[i]=((6./(dx*dx))*(y[i]-2.*y[i+1]+y[i+2])-dp[i-1])/(4.-cp[i-1]);
  }
  M[rows-1]=0.;
  M[rows-2]=dp[rows-3];
  for (i=rows-4;i>-1;i--){
    M[i+1]=dp[i]-cp[i]*M[i+2];
  }
  M[0]=0.;
  for (i=0;i<rows-1;i++)
  {
    coef[0+4*i]=(M[i+1]-M[i])/(6*dx);
    coef[1+4*i]=M[i]/2.;
    coef[2+4*i]=(y[i+1]-y[i])/dx-(M[i+1]+2.*M[i])*dx/6.;
    coef[3+4*i]=y[i];
  }
  //printf("%g\n",coef[rows]);
  //printf("%g\n",*(&coef[0]+rows));
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
