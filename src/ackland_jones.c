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

typedef struct parm_st
{
   unsigned nPairs;
   double rcut;
   int nfiles;
   char* filename;
   char* miscInfo;
} ACKJ_PARMS;

int getRank(int);

typedef struct ackj_store_struct
{
   unsigned ackj_lcs;
   unsigned ackj_cn;
} ACKJ_STORE;

static void ackj_calc(SIMULATE* simulate, ACKJ_PARMS* parms, ACKJ_STORE* ackjData);
static void ackj_write(SIMULATE* simulate, ACKJ_PARMS* parms, ACKJ_STORE* ackjData);

////////////////////////////////////////////////////////////////

ACKJ_PARMS* AcklandJones_parms(ANALYSIS* analysis)
{
   ACKJ_PARMS* parms = ddcMalloc(sizeof(ACKJ_PARMS));
   OBJECT* obj = (OBJECT*) analysis;
   object_get(obj, "nPairs",   &parms->nPairs,   INT,    1, "7");
   object_get(obj, "nfiles",   &parms->nfiles,   INT,    1, "0");
   object_get(obj, "filename", &parms->filename, STRING, 1, "ackjon");
   object_get(obj, "rcut", &parms->rcut, WITH_UNITS, 1, "0.0", "l", NULL);
   double length_convert = units_convert(1.0, NULL, "Angstrom");
   char string[1024];
   sprintf(string, "nPairs = %d; rcut = %f Angstrom;",
           parms->nPairs, parms->rcut * length_convert);
   parms->miscInfo = strdup(string);
   return parms;
}

void AcklandJones_eval(ANALYSIS* analysis)
{
   // empty
}

void AcklandJones_output(ANALYSIS* analysis)
{
   ACKJ_PARMS* parms = analysis->parms;
   SIMULATE* simulate = (SIMULATE*) analysis->parent;
   int nlocal = simulate->system->nlocal;
   unsigned value_blk;
   ACKJ_STORE* ackjData = heapGet(&value_blk);
   heapEndBlock(value_blk, sizeof(ACKJ_STORE) * nlocal);
   ackj_calc(simulate, parms, ackjData);
   ackj_write(simulate, parms, ackjData);
   heapFree(value_blk);
   return;
}

void ackj_calc(SIMULATE* simulate, ACKJ_PARMS* parms, ACKJ_STORE* ackjData)
{
   SYSTEM* sys = simulate->system;
   unsigned nLocal = sys->nlocal;
   unsigned displ_blk;
   RTU_PAIR_FINDER rtupf = rtupf_create(simulate, parms->rcut);
   DISPLACEMENT* displacements = heapGet(&displ_blk);
   unsigned maxNbrs = 0;
   for (unsigned ii = 0; ii < nLocal; ++ii)
   {
      unsigned nNbrs = findNbrs(parms->nPairs, rtupf, ii, displacements);
      assert(nNbrs >= 4 * parms->nPairs);
      maxNbrs = MAX(maxNbrs, nNbrs);
      unsigned N0 = 0, atype = 0;
      unsigned chi_0 = 0, chi_1 = 0, chi_2 = 0, chi_3 = 0;
      unsigned chi_4 = 0, chi_5 = 0, chi_6 = 0, chi_7 = 0;
      double r2_0 = 0., r2_1;
      double cos_theta, delta_bcc, delta_cp, delta_fcc, delta_hcp;

      for (unsigned jj = 0; jj < 6; ++jj)
      {
         DISPLACEMENT* dj = displacements + jj;
         r2_0 += dj->r2;
      }
      r2_0 /= 6.;
      r2_1 = 1.65 * r2_0;

      for (unsigned jj = 0; jj < 4 * parms->nPairs; ++jj)
      {
         DISPLACEMENT* dj = displacements + jj;
         if(dj->r2 < r2_1) N0 = N0 + 1;
      }

      for (unsigned jj = 0; jj < N0; ++jj)
      {
         DISPLACEMENT* dj = displacements + jj;
         for (unsigned kk = jj + 1; kk < N0; ++kk)
         {
            DISPLACEMENT* dk = displacements + kk;
            cos_theta = (dj->x * dk->x + dj->y * dk->y + dj->z * dk->z) / sqrt(dj->r2 * dk->r2);
            if      (-1.001 <= cos_theta && cos_theta < -0.945)chi_0 = chi_0 + 1;
            else if (-0.945 <= cos_theta && cos_theta < -0.915)chi_1 = chi_1 + 1;
            else if (-0.915 <= cos_theta && cos_theta < -0.755)chi_2 = chi_2 + 1;
            else if (-0.755 <= cos_theta && cos_theta < -0.705)chi_3 = chi_3 + 1;
            else if (-0.195 <= cos_theta && cos_theta <  0.195)chi_4 = chi_4 + 1;
            else if ( 0.195 <= cos_theta && cos_theta <  0.245)chi_5 = chi_5 + 1;
            else if ( 0.245 <= cos_theta && cos_theta <  0.795)chi_6 = chi_6 + 1;
            else if ( 0.795 <= cos_theta && cos_theta <= 1.001)chi_7 = chi_7 + 1;
         }
      }

      delta_bcc = 0.35 * chi_4 / (chi_5 + chi_6 - chi_4);
      delta_cp  = fabs((double)chi_6 - 24.) / 24.;
      delta_fcc = 0.61 * (fabs((double)chi_0 + (double)chi_1 - 6.) + (double)chi_2) / 6.;
      delta_hcp = (fabs((double)chi_0 - 3.) +
                   fabs((double)chi_0 + (double)chi_1 + (double)chi_2 + (double)chi_3 - 9.)) / 12.;

      if      (chi_0 == 7)delta_bcc = 0.;
      else if (chi_0 == 6)delta_fcc = 0.;
      else if (chi_0 <= 3)delta_hcp = 0.;

      if      (chi_7 > 0) atype = 0;
      else if (chi_4 < 3.)
      {
         if ((N0 > 13) || (N0 < 11)) atype = 0;
         else atype = 4;
      }
      else if (delta_bcc <= delta_cp)
      {
         if (N0 < 11) atype = 0;
         else atype = 1;
      }
      else if ((N0 > 12) || (N0 < 11)) atype = 0;
      else if (delta_fcc < delta_hcp) atype = 2;
      else atype = 3;

      ackjData[ii].ackj_lcs = atype;
      ackjData[ii].ackj_cn = N0;
   }
   heapEndBlock(displ_blk, sizeof(DISPLACEMENT) * maxNbrs);
   heapFree(displ_blk);
   rtupf_destroy(&rtupf);
}

/** Write out the parameters. */
void ackj_write(SIMULATE* simulate, ACKJ_PARMS* parms, ACKJ_STORE* ackjData)
{
   SYSTEM* sys = simulate->system;
   int nfields = 7;
   char fieldTypes[] = "u u f f f f u";
   char fieldNames[] =
      "checksum label rx ry rz ackj_lcs ackj_cn";
   char fieldUnits[] = "1 1 Angstrom Angstrom Angstrom 1 1" ;
   int lrec = 112;
   CreateSnapshotdir(simulate, NULL);
   char pfilename[1024];
   snprintf(pfilename, 1023, "%s/%s", simulate->snapshotdir, parms->filename);
   PFILE* file = Popen(pfilename, "w", COMM_LOCAL);
   PioSet(file, "recordLength", lrec);
   PioSet(file, "datatype", FIXRECORDASCII);
   PioSet(file, "numberRecords", 0llu);
   PioSet(file, "checksum", CRC32);
   PioSet(file, "nfields", nfields);
   PioSet(file, "field_names", fieldNames);
   PioSet(file, "field_types", fieldTypes);
   PioSet(file, "field_units", fieldUnits);
   PioSet(file, "misc_info", parms->miscInfo);
   if (parms->nfiles > 0) PioSet(file, "ngroup", parms->nfiles);
   if (getRank(0) == 0) write_fileheader(file, simulate, "ackjon" );

   double* rx = sys->collection->state->rx;
   double* ry = sys->collection->state->ry;
   double* rz = sys->collection->state->rz;
   gid_type* label = sys->collection->state->label;
   double length_convert = units_convert(1.0, NULL, "Angstrom");

   for (unsigned jj = 0; jj < sys->nlocal; jj++)
   {
      char line[lrec];
      double x = rx[jj];
      double y = ry[jj];
      double z = rz[jj];
      backInBox_fast(&x, &y, &z);
      x *= length_convert;
      y *= length_convert;
      z *= length_convert;
      unsigned c_lcs = ackjData[jj].ackj_lcs;
      unsigned c_cn  = ackjData[jj].ackj_cn;
      sprintf(line,
              "%08x %12"PRIu64" %14.4f %14.4f %14.4f %8.4f %8.4f",
              0, label[jj], x, y, z, (double)c_lcs, (double)c_cn);
      for (int kk = strlen(line); kk < lrec - 1; ++kk)
         line[kk] = ' ';
      line[lrec-1] = '\n';

      unsigned checksum =
         checksum_crc32_table((unsigned char *)line + 8, lrec - 8);
      char tmp_string[9];
      sprintf(tmp_string, "%08x", checksum);
      memcpy(line, tmp_string, 8);
      Pwrite(line, lrec, 1, file);
   }
   Pclose(file);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
