#include <math.h>
#include "system.h"
#include "collection.h"
#include "state.h"
#include "box.h"
#include "three_algebra.h"
#include "preduce.h"
#include "units.h"
void makeImageParticles(SYSTEM *sys)
{

   double cLen    = units_convert(1.0,NULL, "Ang");
   double rCut = 1.0/cLen; 
   BOX_STRUCT *box = sys->box; 
   STATE *state = sys->collection->state; 
   double* rx = state->rx;
   double* ry = state->ry;
   double* rz = state->rz;
   THREE_MATRIX h = box_get_h(box); 
   THREE_VECTOR *a = (THREE_VECTOR *)&h.xx; 
   THREE_VECTOR *b = (THREE_VECTOR *)&h.yx; 
   THREE_VECTOR *c = (THREE_VECTOR *)&h.zx; 
   THREE_VECTOR aHat = cross(b,c); VSCALE(aHat,1.0/sqrt(VSQ(aHat))); 
   THREE_VECTOR bHat = cross(c,a); VSCALE(bHat,1.0/sqrt(VSQ(bHat))); 
   THREE_VECTOR cHat = cross(a,b); VSCALE(cHat,1.0/sqrt(VSQ(cHat))); 
   THREE_VECTOR d = box_get_diagonal(box);
   THREE_VECTOR corner = box_get_corner(box);
   THREE_VECTOR center = corner; VSVOP(center,+=,0.5,*,d);
   double aL = DOT(aHat,d); 
   double bL = DOT(bHat,d); 
   double cL = DOT(cHat,d); 
   int nion = sys->nion; 
   for (int i=0;i<nion;i++) backInBox(rx+i,ry+i,rz+i);
   int nImage = 0; 
   for (int j=0;j<3;j++)
   {
      THREE_VECTOR hat; 
      double L; 
      THREE_VECTOR *lv; 
      
      if (j==0) {hat = aHat; L = aL; lv= a; }
      if (j==1) {hat = bHat; L = bL; lv= b; }
      if (j==2) {hat = cHat; L = cL; lv= c; }
      int mask = 1; 
      int nT = nion+nImage; 
      for (int i=0;i<nT;i++) 
      {
         THREE_VECTOR diff; 
         diff.x = rx[i]-center.x; 
         diff.y = ry[i]-center.y; 
         diff.z = rz[i]-center.z; 
         double D = DOT(diff,hat); 
         int image  = (D / (0.5*L-rCut))*((box->pbc & mask)>0);
         if (image != 0) 
         {
               resize(nion+nImage+1, 2, state);
               rx[nion+nImage] = rx[i] - image*lv->x;
               ry[nion+nImage] = ry[i] - image*lv->y;
               rz[nion+nImage] = rz[i] - image*lv->z;
               printf("IMAGE %2d %2d %4d %4d %10.4f %10.4f %10.4f\n",j,image,i,nion+nImage,cLen*rx[nion+nImage],cLen*ry[nion+nImage],cLen*rz[nion+nImage]); 
               nImage++; 
         }
       }
     mask *=2; 
   }
}
