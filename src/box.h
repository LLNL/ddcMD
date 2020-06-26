#ifndef BOX_H
#define BOX_H
#include <stdio.h>
#include "three_algebra.h" 
#include "eq.h" 

#ifdef __cplusplus
extern "C" {
#endif

enum BOX_ELEMENTS
{
	BOX_NONE, VOLUME, MINSPAN, HO, HO_PTR, HINV_PTR, CORNER_PTR, HFAC, DHDT,
	VOLUME_FUNCTION_OF_TIME, STRAIN, DEFORMATION_RATE, ROTATION, VOLUME_PTR,
	PBC_PTR, NEAREST_IMAGES, LATTICEVECTORS, RECIP_LATTICEVECTORS, RESET,
	BOX_ENSEMBLE, BOX_TIME, BOX_TIMESTART,
};

enum BOX_CLASS { GENERAL, ORTHORHOMBIC };

typedef struct box_st
{
	char *name;
	char *objclass;
	char *value;
	char *type;
	void *parent; 
	char *ensemble;
	enum BOX_CLASS  itype;
	int pbc;  // periodic boundary condition flag
	int time_dependence; 
	double volume;
	double time; 
	THREE_MATRIX h0, hinv, hfac ;
   THREE_MATRIX dhdt; 
	EQTARGET  *Veq; 
	struct {EQTARGET *xx,*xy,*xz,*yx,*yy,*yz,*zx,*zy,*zz;} dhdt_hFunc;
	THREE_MATRIX deformationRate;
	THREE_MATRIX rotationMatrix;
	THREE_VECTOR reducedcorner;
	THREE_VECTOR corner;
	THREE_VECTOR bodyDiag; // longest body diagonal
	double r2Inner; // r^2 of inscribed sphere
	double r2Outer; // r^2 of circumsphere
	int lvSize;     // number of elements in lv
	int lvCapacity; // max capacity of lv
	THREE_INT* lv;   // lattice vectors to image cells that overlap Wigner-Seitz cell.
	unsigned nAffineUpdates;
	THREE_MATRIX* affineUpdates;
	
} BOX_STRUCT;

BOX_STRUCT *box_init(void *parent, char *name);
void box_read(BOX_STRUCT *box,int mode,FILE *file);
void box_read(BOX_STRUCT *box,int mode,FILE *file);
void box_write(BOX_STRUCT*box, FILE*file);
void box_put(BOX_STRUCT*box, int cmd, void *ptr);
BOX_STRUCT *box_getBox(BOX_STRUCT *box) ;
int box_getType(BOX_STRUCT *box) ;
THREE_VECTOR box_get_diagonal(BOX_STRUCT *box) ;
THREE_MATRIX box_get_h(BOX_STRUCT *box) ;
THREE_MATRIX box_get_dhdt(BOX_STRUCT *box) ;
THREE_VECTOR box_get_corner(BOX_STRUCT *box) ;
THREE_VECTOR box_get_reducedcorner(BOX_STRUCT *box) ;
THREE_VECTOR box_get_urcorner(BOX_STRUCT *box) ;
int  box_get_boundary_conditions(BOX_STRUCT *box) ;
double  box_get_volume(BOX_STRUCT *box) ;
double box_get_minspan(BOX_STRUCT  *box);
THREE_VECTOR box_get_boundingbox(BOX_STRUCT *box) ;
void box_get(BOX_STRUCT*box, int cmd, void *ptr);
unsigned box_newAffineUpdateIndex(BOX_STRUCT* box);
THREE_MATRIX box_getAffineUpdateHfac(BOX_STRUCT* box, unsigned index);
void boxPrescriptiveTimeParse(BOX_STRUCT *box);
THREE_MATRIX boxPrescriptiveTime(BOX_STRUCT *box,double newTime);

#ifdef __cplusplus
}
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
