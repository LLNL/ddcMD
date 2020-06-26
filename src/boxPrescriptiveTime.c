#include "box.h"
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "ddcMalloc.h"
#include "gid.h"
#include "eq.h"
#include "system.h"

void boxPrescriptiveTimeParse(BOX_STRUCT *box)
{
   box->time_dependence=BOX_NONE; 
   char **u = NULL; 
   EQTARGET *u_eq[9]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}; 
   int nelements = object_getv((OBJECT *) box, "dudt", (void *)&u, STRING,IGNORE_IF_NOT_FOUND);
   if (u!=NULL) 
   {
      box->time_dependence=STRAIN; 
      box->dhdt_hFunc.xx = box->dhdt_hFunc.yx = box->dhdt_hFunc.zx = NULL;
      box->dhdt_hFunc.xy = box->dhdt_hFunc.yy = box->dhdt_hFunc.zy = NULL;
      box->dhdt_hFunc.xz = box->dhdt_hFunc.yz = box->dhdt_hFunc.zz = NULL;
      if (nelements == 0) 
      {
         u_eq[6] = u_eq[3] = u_eq[0]=eq_parse("0.0","1/t","t");
         u_eq[7] = u_eq[4] = u_eq[1]=eq_parse("0.0","1/t","t");
         u_eq[8] = u_eq[5] = u_eq[2]=eq_parse("0.0","1/t","t");
      }
      if (nelements == 1) 
      {
         u_eq[6] = u_eq[3] = u_eq[0]=eq_parse(u[0],"1/t","t");
         u_eq[7] = u_eq[4] = u_eq[1]=eq_parse(u[0],"1/t","t");
         u_eq[8] = u_eq[5] = u_eq[2]=eq_parse(u[0],"1/t","t");
      }
      if (nelements == 2) 
      {
         u_eq[6] = u_eq[3] = u_eq[0]=eq_parse(u[0],"1/t","t");
         u_eq[7] = u_eq[4] = u_eq[1]=eq_parse(u[1],"1/t","t");
         u_eq[8] = u_eq[5] = u_eq[2]=eq_parse(u[1],"1/t","t");
      }
      if (nelements == 3) 
      {
         u_eq[6] = u_eq[3] = u_eq[0]=eq_parse(u[0],"1/t","t");
         u_eq[7] = u_eq[4] = u_eq[1]=eq_parse(u[1],"1/t","t");
         u_eq[8] = u_eq[5] = u_eq[2]=eq_parse(u[2],"1/t","t");
      }
      if (nelements == 9) 
      {
         u_eq[0]=eq_parse(u[0],"1/t","t");
         u_eq[1]=eq_parse(u[1],"1/t","t");
         u_eq[2]=eq_parse(u[2],"1/t","t");
         u_eq[3]=eq_parse(u[3],"1/t","t");
         u_eq[4]=eq_parse(u[4],"1/t","t");
         u_eq[5]=eq_parse(u[5],"1/t","t");
         u_eq[6]=eq_parse(u[6],"1/t","t");
         u_eq[7]=eq_parse(u[7],"1/t","t");
         u_eq[8]=eq_parse(u[8],"1/t","t");
      }
      box->dhdt_hFunc.xx =  u_eq[0]; 
      box->dhdt_hFunc.xy =  u_eq[1]; 
      box->dhdt_hFunc.xz =  u_eq[2]; 
      box->dhdt_hFunc.yx =  u_eq[3]; 
      box->dhdt_hFunc.yy =  u_eq[4]; 
      box->dhdt_hFunc.yz =  u_eq[5]; 
      box->dhdt_hFunc.zx =  u_eq[6]; 
      box->dhdt_hFunc.zy =  u_eq[7]; 
      box->dhdt_hFunc.zz =  u_eq[8]; 
      for (int ii=0; ii<nelements; ++ii) ddcFree(u[ii]);
      ddcFree(u);
      u = NULL;
     return; 
   }

   char *string; 
   nelements = object_get((OBJECT *) box, "Veq", &string, LITERAL,1,"");
   if (*string != '\0') 
   {
      prune_spaces(string); 
      box->time_dependence=VOLUME_FUNCTION_OF_TIME; 
      box->Veq = eq_parse(string,"l^3","t");
      if (box->Veq != NULL && box->Veq->function) (box->Veq->function) (box->time, (void *)box->Veq);
      ddcFree(string);
      return; 
   }

   nelements = object_get(
         (OBJECT*)box, "deformationRate", &(box->deformationRate),
         WITH_UNITS, 9,	"0 0 0  0 0 0  0 0 0", "1/t", NULL);
   if (! matrix_equal(box->deformationRate, mzero) ) box->time_dependence = DEFORMATION_RATE;

   nelements = object_get(
         (OBJECT*)box, "rotationMatrix", &(box->rotationMatrix),
         DOUBLE, 9, "0 0 0  0 0 0  0 0 0");
   if (! matrix_equal(box->rotationMatrix, mzero) ) box->time_dependence = ROTATION;

}
THREE_MATRIX boxPrescriptiveTime(BOX_STRUCT *box,double newTime)
{
  double time = box->time; 
  THREE_MATRIX h = box->h0; 
   switch(box->time_dependence)
   {
      case STRAIN:
         h.xx *= exp(box->dhdt_hFunc.xx->integral(time,newTime,(void*)box->dhdt_hFunc.xx)) ; 
         h.yx *= exp(box->dhdt_hFunc.yx->integral(time,newTime,(void*)box->dhdt_hFunc.yx)) ; 
         h.zx *= exp(box->dhdt_hFunc.zx->integral(time,newTime,(void*)box->dhdt_hFunc.zx)) ; 
         h.xy *= exp(box->dhdt_hFunc.xy->integral(time,newTime,(void*)box->dhdt_hFunc.xy)) ; 
         h.yy *= exp(box->dhdt_hFunc.yy->integral(time,newTime,(void*)box->dhdt_hFunc.yy)) ; 
         h.zy *= exp(box->dhdt_hFunc.zy->integral(time,newTime,(void*)box->dhdt_hFunc.zy)) ; 
         h.xz *= exp(box->dhdt_hFunc.xz->integral(time,newTime,(void*)box->dhdt_hFunc.xz)) ; 
         h.yz *= exp(box->dhdt_hFunc.yz->integral(time,newTime,(void*)box->dhdt_hFunc.yz)) ; 
         h.zz *= exp(box->dhdt_hFunc.zz->integral(time,newTime,(void*)box->dhdt_hFunc.zz)) ; 
         break; 
      case VOLUME_FUNCTION_OF_TIME:
         {
            gid_type nglobal=system_getNglobal(NULL);
            double a = cbrt(nglobal*box->Veq->function(newTime,box->Veq)/box->volume);
            MATSCALE(h,a); 
         }
         break; 
      case DEFORMATION_RATE:
         { //limit scope
            double dt = newTime-time;
            THREE_MATRIX sum = I_3x3;
            THREE_MATRIX xx = matrix_sadd(dt, box->deformationRate, mzero);
            THREE_MATRIX term = xx;
            unsigned kk = 0;
            double invFactorial = 1;
            do
            {
               ++kk;
               invFactorial /= (1.0*kk);
               sum = matrix_sadd(invFactorial, term, sum);
               term = matrix_matrix(xx, term);
            }
            while ( ! matrix_equal_tol(term, mzero, 1e-15) );

            h = matrix_matrix(h, sum);
         } // limit scope
         break;
      case ROTATION:
         h = matrix_matrix(box->rotationMatrix, h);
         break;
   }
   return h; 
}
