
////////////////////////////////////////////////////////////////////////
//
//  This code was re-written in r1161 to match the implementation in Rob
//  Rudd's MD code.  This implementation substantially reduces the
//  effect of thermal noise.  It also offers an improved parameter set.
//  See the sourceforge wiki or the svn logs for more information.
//
////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "preduce.h"
#include "object.h"
#include "pio.h"
#include "io.h"
#include "simulate.h"
#include "system.h"
#include "heap.h"
#include "ddcMalloc.h"
#include "crc32.h"
#include "units.h"

#include "rtuPairFinder.h"

int getRank(int);

typedef struct parm_st
{
   unsigned nPairs;
   double rcut;
   double V0;
   double cen_th;
   int nfiles;
   char* filename;
   char* miscInfo;
} CSYM_PARMS;

typedef struct csym_store_struct
{
   double csym_cm;
   double scale;
} CSYM_STORE;

static void csym_calc(SIMULATE* simulate, CSYM_PARMS* parms,
		      CSYM_STORE* csymData);
static void csym_write(SIMULATE* simulate, CSYM_PARMS* parms,
		       CSYM_STORE* csymData);

////////////////////////////////////////////////////////////////


CSYM_PARMS* centroSym_parms(ANALYSIS* analysis)
{
   CSYM_PARMS* parms = ddcMalloc(sizeof(CSYM_PARMS));
   OBJECT* obj = (OBJECT*) analysis;
   object_get(obj, "nPairs",   &parms->nPairs,   INT,    1, "7");
   object_get(obj, "nfiles",   &parms->nfiles,   INT,    1, "0");
   object_get(obj, "filename", &parms->filename, STRING, 1, "censym");
   object_get(obj, "cen_th",   &parms->cen_th,   WITH_UNITS, 1, "0.0", "l^2", NULL);
   object_get(obj, "rcut",     &parms->rcut,     WITH_UNITS, 1, "0.0", "l",   NULL);
   object_get(obj, "V0",       &parms->V0,       WITH_UNITS, 1, "0.0", "l^3", NULL);

   double length_convert = units_convert(1.0, NULL, "Angstrom"); 
   char string[1024];
   sprintf(string,
	   "nPairs = %d; cen_th = %f Angstrom^2; rcut = %f Angstrom; V0 = %f Angstrom^3;",
	   parms->nPairs, parms->cen_th*length_convert*length_convert,
	   parms->rcut*length_convert, parms->V0*pow(length_convert,3));
   
   parms->miscInfo = strdup(string);
   return parms;
}

void centroSym_eval(ANALYSIS* analysis)
{
   // empty
}

/** Evaluate the instantaneous centrosymmetry order parameter at all
 *  atomic sites in the system and write the results to a file.  */
void centroSym_output(ANALYSIS* analysis)
{
   CSYM_PARMS* parms = analysis->parms;
   SIMULATE* simulate = (SIMULATE*) analysis->parent; 
   int nlocal = simulate->system->nlocal; 

   unsigned value_blk;
   CSYM_STORE* csymData = heapGet(&value_blk);
   heapEndBlock(value_blk, sizeof(CSYM_STORE) * nlocal);

   csym_calc(simulate, parms, csymData); 
   csym_write(simulate, parms, csymData); 
   heapFree(value_blk);

   return;
}


/**
 * Compute local centrosymmetry order parameter for nlocal atoms.  The
 * algorithm identifies PAIRS j,k that are both neighbors around the
 * central atom and uses the deviation of the jk center of mass from an
 * "origin" to quantify centrosymmetric order.
 *                                                                   
 * References: C.L. Kelchner, S.J. Plimpton, J.C. Hamilton,          
 * "Dislocation nucleation and defect structure during surface       
 * indentation," Phys. Rev. B58, 11085-11088 (1998),                 
 * (first use of a centrosymmetry parameter)                         
 *                                                                   
 * and following http://mt.seas.upenn.edu/Archive/\                  
 *                    Graphics/A/Doc/CentralSymmetry.pdf,            
 * "Central Symmetry Parameter" by Ju Li, August 30, 2003.           
 *
 * We also include enhancements due to R.E. Rudd that mitigate the
 * effect of thermal noise by using the center of mass of the central
 * atom and its neighbors as the origin instead of the actual atom
 * position.
 *
 *                                                                   
 * This version of centrosymmetry writes output files every outputrate
 * steps.  The eval function is empty.  There is no capability to
 * accumulate and average centrosymmetry data.
 *
 *
 * This function assumes that findNbrs always returns n4 =
 * 4*parms->nPairs neighbor atoms.  In the event that an atom doesn't
 * have enough neighbors in range, synthetic data is used with zero
 * displacement, but a large distance so that it will sort to the end of
 * the displacement list.
 *
 * Once the neighbors are found we calculate the center of mass (not
 * really the center of mass because we don't account for possible mass
 * differences) of the central atom and its nearest n2 = 2*parms->nPairs
 * atoms.  Note that the center of mass position is expressed in
 * corrdinates relative to the position of the central atom.  We then
 * loop over the neighbor atoms, starting with the closest, and find
 * pairs that minimize the centrosymmetry deviation with respect to
 * the center of mass:
 *
 *   c = | dj-rcm + dk-rcm |^2
 *
 * When we have found enough pairs we stop.  Since we have n4 neighbors
 * to choose from, we can always find enough pairs even though some may
 * involve pairs with synthetic atoms that sit at the center of mass.
 *
 *  One other detail: The displacements scratch array is obtained from
 *  the scratch heap, but we don't close the block until the end of the
 *  routine.  Since the scratch heap is large compared to the number of
 *  nbrs we should find is practically impossible to overrun the array.
 *  (And the user can make the heap bigger at run time if needed.)
 *  We still keep track of how many nbrs we found so that the heap code
 *  can be sure that we didn't overrun when heapEndBlock is called.
 */
void csym_calc(SIMULATE* simulate, CSYM_PARMS* parms, CSYM_STORE* csymData)
{
   SYSTEM* sys = simulate->system;
   unsigned nLocal = sys->nlocal;

   unsigned displ_blk;
   RTU_PAIR_FINDER rtupf = rtupf_create(simulate, parms->rcut);
   DISPLACEMENT* displacements = heapGet(&displ_blk); 

   const unsigned n1 = parms->nPairs;
   const unsigned n2 = 2 * parms->nPairs;
   const unsigned n4 = 4 * parms->nPairs;

   unsigned maxNbrs = 0;
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      unsigned nNbrs = findNbrs(parms->nPairs, rtupf, ii, displacements);
      assert(nNbrs >= n4);
      maxNbrs = MAX(maxNbrs, nNbrs);  // track amount of heap used.

      // calculate rcm relative to atom ii.
      THREE_VECTOR rcm = vzero;
      for (unsigned jj=0; jj<n2; ++jj)
      {
	 DISPLACEMENT* dj = displacements + jj;
	 rcm.x += dj->x;
	 rcm.y += dj->y;
	 rcm.z += dj->z;
      }
      rcm.x /= (n2 + 1.0);
      rcm.y /= (n2 + 1.0);
      rcm.z /= (n2 + 1.0);
      
      unsigned bprocessed[n4];
      for (unsigned jj=0; jj<n4; ++jj) 
	 bprocessed[jj]=0;	 

      // Loop over neighbors, j, starting with the closest and working
      // outwards.  Pair the "closest" with its "best" centrosymmetric
      // partner, k>j, and sum their centrosymmetry deviation
      // parameters.  Mark the pair as having been counted/processed.
      // Keep going until enough pairs have been identified.
      unsigned nPairs = 0;
      double u2_sum = 0.0;
      double scale = 0.0;
      for (unsigned jj=0; jj<n2; ++jj)
      {
	 if (nPairs == n1) break;
	 if (bprocessed[jj] == 1)
	    continue;  // skip processed atoms
      
	 bprocessed[jj]=1; 
	 DISPLACEMENT* dj = displacements+jj; 
	 unsigned kMin = jj;
	 double u2Min = 4.0 * rtupf.rcut * rtupf.rcut;
	 for (unsigned kk=jj+1; kk<n4; ++kk)
	 {
	    if (bprocessed[kk] != 0)
	       continue;
	    
	    DISPLACEMENT* dk = displacements+kk; 
	    
	    THREE_VECTOR u;
	    u.x = dj->x + dk->x - 2.0*rcm.x;
	    u.y = dj->y + dk->y - 2.0*rcm.y; 
	    u.z = dj->z + dk->z - 2.0*rcm.z;
	    double u2 = u.x*u.x + u.y*u.y + u.z*u.z;		
	    if (u2 < u2Min) 
	    {
	       u2Min = u2;
	       kMin = kk;
	    }
	 }
	 assert(kMin != jj);
	 bprocessed[kMin] = 1; 
	 nPairs += 1;
	 u2_sum += u2Min;
	 DISPLACEMENT* dk = displacements+kMin; 
	 // can't get scale from dj->r2 since that distance might be
	 // artificial.
	 scale += dj->x*dj->x + dj->y*dj->y + dj->z*dj->z;
	 scale += dk->x*dk->x + dk->y*dk->y + dk->z*dk->z;
      }

      csymData[ii].csym_cm = u2_sum; 
      csymData[ii].scale = scale;
   }

   heapEndBlock(displ_blk, sizeof(DISPLACEMENT) * maxNbrs);
   heapFree(displ_blk);
   rtupf_destroy(&rtupf);
}


/** Write out the order parameters. */
void csym_write(SIMULATE* simulate, CSYM_PARMS* parms, CSYM_STORE* csymData)
{
   SYSTEM* sys = simulate->system;
   
   int nfields = 7;
   char fieldTypes[] = "u u f f f f f";
   char fieldNames[] =
      "checksum label rx ry rz csym_cm csym_scaled";
   char fieldUnits[] = "1 1 Angstrom Angstrom Angstrom Angstrom^2 1 " ;
   int lrec = 112; // 105 bytes per record and some padding.

   CreateSnapshotdir(simulate, NULL);
   char pfilename[1024];
   snprintf(pfilename, 1023, "%s/%s",simulate->snapshotdir, parms->filename);
   PFILE* file = Popen(pfilename, "w", COMM_LOCAL);
   PioSet(file, "recordLength", lrec);
   PioSet(file, "datatype", FIXRECORDASCII);
   // we don't know how many records there will be.
   PioSet(file, "numberRecords", 0llu);
   PioSet(file, "checksum", CRC32);
   PioSet(file, "nfields", nfields);
   PioSet(file, "field_names", fieldNames);
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "field_units", fieldUnits);
   PioSet(file, "misc_info", parms->miscInfo);
   if (parms->nfiles > 0) PioSet(file, "ngroup", parms->nfiles);
   if (getRank(0) == 0) write_fileheader(file, simulate, "censym" );
   
   double* rx = sys->collection->state->rx;
   double* ry = sys->collection->state->ry;
   double* rz = sys->collection->state->rz;
   gid_type* label = sys->collection->state->label;

   double length_convert = units_convert(1.0, NULL, "Angstrom"); 
   double length_convert2 = length_convert * length_convert;

   for (unsigned jj = 0; jj < sys->nlocal; jj++)
   {
      if (csymData[jj].csym_cm >= parms->cen_th)
      {
	 char line[lrec];
	 double x = rx[jj];
	 double y = ry[jj];
	 double z = rz[jj];
	 backInBox_fast(&x, &y, &z);
	 x *= length_convert;
	 y *= length_convert;
	 z *= length_convert;
	 double c_cm = csymData[jj].csym_cm * length_convert2;
	 double c_s = csymData[jj].csym_cm/(2 * csymData[jj].scale);
	 sprintf(line,
                 "%08x %12"PRIu64" %14.4f %14.4f %14.4f %18.8f %18.8e",
		 0, label[jj], x, y, z, c_cm, c_s);
	 for (int kk=strlen(line); kk<lrec-1; ++kk)
	    line[kk] = ' ';
	 line[lrec-1] = '\n';

	 unsigned checksum =
                 checksum_crc32_table((unsigned char *)line+8,lrec-8);
	 char tmp_string[9];
	 sprintf(tmp_string,"%08x",checksum);
	 memcpy(line, tmp_string, 8);
	 Pwrite(line, lrec, 1, file);
      } 
   }
   Pclose(file); 
}
