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
   double rfcut;
   unsigned NNs;
   int nfiles;
   char* filename;
   char* miscInfo;
} QUATERNION_PARMS;

typedef struct quaternion_store_struct
{
   double quaternion_0;
   double quaternion_1;
   double quaternion_2;
   double quaternion_3;
   double quaternion_h;
} QUATERNION_STORE;

static void quaternion_calc(SIMULATE* simulate, QUATERNION_PARMS* parms,
                            QUATERNION_STORE* quaternionData);
static void quaternion_write(SIMULATE* simulate, QUATERNION_PARMS* parms,
                             QUATERNION_STORE* quaternionData);

////////////////////////////////////////////////////////////////

QUATERNION_PARMS* quaternion_parms(ANALYSIS* analysis)
{
   QUATERNION_PARMS* parms = ddcMalloc(sizeof(QUATERNION_PARMS));
   OBJECT* obj = (OBJECT*) analysis;
   object_get(obj, "nPairs",   &parms->nPairs,   INT,    1, "7");
   object_get(obj, "rcut", &parms->rcut, WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "rfcut", &parms->rfcut, DOUBLE, 1, "1.65");
   object_get(obj, "NNs", &parms->NNs, INT, 1, "8");
   object_get(obj, "nfiles",   &parms->nfiles,   INT,    1, "0");
   object_get(obj, "filename", &parms->filename, STRING, 1, "quaternion");
   double length_convert = units_convert(1.0, NULL, "Angstrom");
   char string[1024];
   sprintf(string, "nPairs = %d; rcut = %f Angstrom; rfcut = %f;",
           parms->nPairs, parms->rcut * length_convert, parms->rfcut);
   parms->miscInfo = strdup(string);
   return parms;
}

void quaternion_eval(ANALYSIS* analysis)
{
   // empty
}

void quaternion_output(ANALYSIS* analysis)
{
   QUATERNION_PARMS* parms = analysis->parms;
   SIMULATE* simulate = (SIMULATE*) analysis->parent;
   int nlocal = simulate->system->nlocal;
   unsigned value_blk;
   QUATERNION_STORE* quaternionData = heapGet(&value_blk);
   heapEndBlock(value_blk, sizeof(QUATERNION_STORE) * nlocal);
   quaternion_calc(simulate, parms, quaternionData);
   quaternion_write(simulate, parms, quaternionData);
   heapFree(value_blk);
   return;
}

void quaternion_calc(SIMULATE* simulate, QUATERNION_PARMS* parms, QUATERNION_STORE* quaternionData)
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
      double r2_0 = 0., r2_1;
      double dx, dy, dz, dr = 0., djkn = 0., djkm = 0.;
      double cos_theta, theta, phi, psi;
      int jn = -1, jm = -1, kn = -1, km = -1;
      unsigned N0 = 0, NNs = 0;
      double nxyz = -3., mxyz = -3., jxyz = -3., kxyz = -3.;
      double nx, ny, nz, mx, my, mz, px, py, pz, pn = 0., qx, qy, qz, qn = 0.;
      //double Q0 = -0.1; 
      //double QH = -0.1; 
      double QR = -0.1, QG = -0.1, QB = -0.1,  factor = 1. - 1e-5;

      for (unsigned jj = 0; jj < 6; ++jj)
      {
         DISPLACEMENT* dj = displacements + jj;
         r2_0 += dj->r2;
      }
      r2_0 /= 6.;
      r2_1 = parms->rfcut * r2_0;

      for (unsigned jj = 0; jj < 4 * parms->nPairs; ++jj)
      {
         DISPLACEMENT* dj = displacements + jj;
         if(dj->r2 < r2_1) N0 = N0 + 1;
      }

      for (unsigned jj = 0; jj < N0; ++jj)
      {
         DISPLACEMENT* djj = displacements + jj;
         for (unsigned kk = 0; kk < N0; ++kk)
         {
            DISPLACEMENT* dkk = displacements + kk;
            cos_theta = (djj->x * dkk->x + djj->y * dkk->y + djj->z * dkk->z) / sqrt(djj->r2 * dkk->r2);
            if (-1.001 <= cos_theta && cos_theta < -0.945)
            {
               NNs = NNs + 1;
               dr = sqrt((djj->x - dkk->x) * (djj->x - dkk->x) +
                         (djj->y - dkk->y) * (djj->y - dkk->y) +
                         (djj->z - dkk->z) * (djj->z - dkk->z));
               dx = (djj->x - dkk->x) / dr;
               dy = (djj->y - dkk->y) / dr;
               dz = (djj->z - dkk->z) / dr;
               if ((( dx + dy + dz) >= nxyz))
               {
                  nxyz =  dx + dy + dz;
                  jn = jj;
                  kn = kk;
               }
               if (((-dx + dy + dz) >= mxyz))
               {
                  mxyz = -dx + dy + dz;
                  jm = jj;
                  km = kk;
               }
               if ((( dx - dy + dz) >= jxyz))
               {
                  jxyz =  dx - dy + dz;
               }
               if ((( dx + dy - dz) >= kxyz))
               {
                  kxyz =  dx + dy - dz;
               }
            }
         }
      }

      //if (NNs == 8 || NNs == 12)
      if (NNs == parms->NNs)
      {
         DISPLACEMENT* djn = displacements + jn;
         DISPLACEMENT* djm = displacements + jm;
         DISPLACEMENT* dkn = displacements + kn;
         DISPLACEMENT* dkm = displacements + km;
         djkn = sqrt((djn->x - dkn->x) * (djn->x - dkn->x) +
                     (djn->y - dkn->y) * (djn->y - dkn->y) +
                     (djn->z - dkn->z) * (djn->z - dkn->z));
         djkm = sqrt((djm->x - dkm->x) * (djm->x - dkm->x) +
                     (djm->y - dkm->y) * (djm->y - dkm->y) +
                     (djm->z - dkm->z) * (djm->z - dkm->z));
         nx = (djn->x - dkn->x) / djkn;
         ny = (djn->y - dkn->y) / djkn;
         nz = (djn->z - dkn->z) / djkn;
         mx = (djm->x - dkm->x) / djkm;
         my = (djm->y - dkm->y) / djkm;
         mz = (djm->z - dkm->z) / djkm;
         px = ny * mz - nz * my;
         py = nz * mx - nx * mz;
         pz = nx * my - ny * mx;
         pn = sqrt(px * px + py * py + pz * pz);
         px = px / pn;
         py = py / pn;
         pz = pz / pn;
         qx = ny * pz - nz * py;
         qy = nz * px - nx * pz;
         qz = nx * py - ny * px;
         qn = sqrt(qx * qx + qy * qy + qz * qz);
         qx = qx / qn;
         qy = qy / qn;
         qz = qz / qn;

         theta = acos(factor * (nx + ny + nz) / sqrt(3.));
         if (theta == 0.)
         {
            phi = 0.;
            psi = acos(factor * (0. - py + pz) / sqrt(2.));
         }
         else
         {
            phi = asin(factor * (0. - ny + nz) / (sqrt(2.) * sin(theta)));
            psi = asin(factor * (px + py + pz) / (sqrt(3.) * sin(theta)));
         }

         //Q0 =       cos(theta / 2.) * cos((phi + psi) / 2.); //nx * mx + ny * my + nz * mz;
         QR = (1. + sin(theta / 2.) * cos((phi - psi) / 2.)) / 2.; //asin(-pz);
         QG = (1. + sin(theta / 2.) * sin((phi - psi) / 2.)) / 2.; //atan2(nz, qz);
         QB = (1. + cos(theta / 2.) * sin((phi + psi) / 2.)) / 2.; //atan2(px, py);
         //QR =     sin(theta/2.)*cos((phi-psi)/2.);
         //QG =     sin(theta/2.)*sin((phi-psi)/2.);
         //QB =     cos(theta/2.)*sin((phi+psi)/2.);
         //QR = (theta+pi/2.)/pi;
         //QG = (phi+pi/2.)/pi;
         //QB = (psi+pi/2.)/pi;

         //QH = 0.;
         //if (QR >  QG && QG >  QB)QH =      (QG - QB) / (QR - QB);
         //if (QG >= QR && QR >  QB)QH = 2. - (QR - QB) / (QG - QB);
         //if (QG >  QB && QB >= QR)QH = 2. + (QB - QR) / (QG - QR);
         //if (QB >= QG && QG >  QR)QH = 4. - (QG - QR) / (QB - QR);
         //if (QB >  QR && QR >= QG)QH = 4. + (QR - QG) / (QB - QG);
         //if (QR >= QB && QB >= QG)QH = 6. - (QB - QG) / (QR - QG);
      }

      quaternionData[ii].quaternion_0 = (QR + QG + QB) / 3.;
      quaternionData[ii].quaternion_1 = QR;
      quaternionData[ii].quaternion_2 = QG;
      quaternionData[ii].quaternion_3 = QB;
      quaternionData[ii].quaternion_h = (QR * QG * QB); //QH / 6.;
   }

   heapEndBlock(displ_blk, sizeof(DISPLACEMENT) * maxNbrs);
   heapFree(displ_blk);
   rtupf_destroy(&rtupf);
}

/** Write out the parameters. */
void quaternion_write(SIMULATE* simulate, QUATERNION_PARMS* parms, QUATERNION_STORE* quaternionData)
{
   SYSTEM* sys = simulate->system;
   int nfields = 10;
   char fieldTypes[] = "u u f f f f f f f f";
   char fieldNames[] =
      "checksum label rx ry rz quaternion_0 quaternion_1 quaternion_2 quaternion_3 quaternion_h";
   char fieldUnits[] = "1 1 Angstrom Angstrom Angstrom 1 1 1 1 1" ;
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
   if (getRank(0) == 0) write_fileheader(file, simulate, "quaternion" );

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
      double c_0 = quaternionData[jj].quaternion_0;
      double c_1 = quaternionData[jj].quaternion_1;
      double c_2 = quaternionData[jj].quaternion_2;
      double c_3 = quaternionData[jj].quaternion_3;
      double c_h = quaternionData[jj].quaternion_h;
      sprintf(line,
              "%08x %12"PRIu64" %14.4f %14.4f %14.4f %8.4f %8.4f %8.4f %8.4f %8.4f",
              0, label[jj], x, y, z, c_0, c_1, c_2, c_3, c_h);
      for (int kk = strlen(line); kk < lrec - 1; ++kk)
         line[kk] = ' ';
      line[lrec - 1] = '\n';

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
