#include "ptiming.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "pio.h"
#include "hardwareInfo.h"
#include "mpiUtils.h"
#include "utilities.h"
#include "rectimer.h"

static PTIMING timing[NTIMING];
static char threadScheduleMethodName[30];

PTIMING *profile_get_timing(void)
{
	return timing; 
}

double profileN(int n)
{
	static double Nsum=0.0;
	static int Ncall=0.0;

	if (n == -2)
	{
		Nsum = Ncall = 0;
		return 0.0;
	}
	if (n == -1) return (Nsum/Ncall);
	Nsum += n;
	Ncall++;
	return Nsum/Ncall;
}

double profileN2Search(int n2)
{
   static double Nsum2=0.0;
   static int Ncall=0.0;

   if (n2 == -2)
   {
      Nsum2 = Ncall = 0;
      return 0.0;
   }
   if (n2 == -1) { return (Nsum2/Ncall);}
   Nsum2 += n2;
   Ncall++;
   return Nsum2/Ncall;
}
double profileN2(int n2)
{
   static double Nsum2=0.0;
   static int Ncall=0.0;

   if (n2 == -2)
   {
      Nsum2 = Ncall = 0;
      return 0.0;
   }
   if (n2 == -1) { return (Nsum2/Ncall);}
   Nsum2 += n2;
   Ncall++;
   return Nsum2/Ncall;
}
double profileN2R(int n2)
{
   static double Nsum2=0.0;
   static int Ncall=0.0;

   if (n2 == -2)
   {
      Nsum2 = Ncall = 0;
      return 0.0;
   }
   if (n2 == -1) { return (Nsum2/Ncall);}
   Nsum2 += n2;
   Ncall++;
   return Nsum2/Ncall;
}

double profileNBox(int nbox)
{
   static double NBoxSum=0.0;
   static int Ncall=0.0;

   if (nbox == -2)
   {
      NBoxSum = Ncall = 0;
      return 0.0;
   }
   if (nbox == -1) 
   {
      if (Ncall == 0) return 0.0;
      return (NBoxSum/Ncall);
   }
   NBoxSum += nbox;
   Ncall++;
   return NBoxSum/Ncall;
}

char* profileThreadScheduleMethodName(char* st)
{
   if (st == NULL)
      return threadScheduleMethodName;
   else
      strcpy(threadScheduleMethodName,st);
   return NULL;
}
//identify thread memory

double profileThrdImbalance(double m)
{
   static double MTsum;
   static int Ncall;

   if (m == -2.0)
   {
      MTsum = Ncall = 0;
      return 0.0;
   }
   if (m == -1.0) return (MTsum/Ncall);
   MTsum += m;
   Ncall++;
   return MTsum/Ncall;
}

double profileMT(int m)
{
   static int numSchedule;
   static int Ncall;
   static double memNeed;

   if (m == -2)
   {
      numSchedule = Ncall = memNeed = 0;
      return 0.0;
   }
   if (m == -1) return (memNeed/Ncall);
   if (m == -3) return (numSchedule);
   if (m == -4) return (numSchedule++);
   memNeed += m;
   Ncall++;
   return memNeed/Ncall;
}

void profileSetLoopTime(double time)
{
	timing[LOOP].elapsetotal=time;
}

double profileGetLoopTime(void)
{
	return timing[LOOP].elapsetime; 
}

// get time associated to local operations
// and contributing to parallel load balance/inbalance
// Note: Use .elapsetotal as it is not reset to 0 at START of timer
// and timers may be running
double profileGetLocalLoadElapsetime(int flag)
{
   static int ncalls=0;
   static double previous_elapsetotal=0.;
   static double previous_elapsetotalc[5]={0.,0.,0.,0.,0.};
   
   static double previous_ncv=0;
   static double previous_ncu=0;
   static double previous_ncf=0;
   
   double time = 0.;
   double deltat = 0.;

   switch( flag ){
      case 0:
         {
            deltat    = timing[MDSTEP].elapsetotal - previous_elapsetotal;
            previous_elapsetotal += deltat;

            double timec[5] = {timing[BARRIER1].elapsetotal,
                               timing[BARRIER2].elapsetotal,
                               timing[WAIT_FPAIRS].elapsetotal,
                               timing[EVAL_ETYPE].elapsetotal,
                               timing[COMMCOMM].elapsetotal};

            for(short i=0;i<5;i++){
               timec[i]                 -= previous_elapsetotalc[i];
               previous_elapsetotalc[i] += timec[i];
            }

            time = deltat;
            for(short i=0;i<5;i++)time-=timec[i];
            
            assert( time>=0. );
         }
         break;
      
      case 1:
         time = timing[P_FORCE].elapsetotal - previous_elapsetotal;
         previous_elapsetotal += time;
         break;

      case 2:
         {      
            double ncf = timing[P_FORCE].ncalls;
            double ncv = timing[VERLET].ncalls;
            double ncu = timing[UPDATEPAIR].ncalls;
            
            double dncv = ncv - previous_ncv;
            double dncu = ncu - previous_ncu;
            double dncf = ncf - previous_ncf;
            
            double timec[3] = {timing[P_FORCE].elapsetotal,
                               timing[UPDATEPAIR].elapsetotal,
                               timing[VERLET].elapsetotal};

            for(short i=0;i<3;i++){
               timec[i]                 -= previous_elapsetotalc[i];
               previous_elapsetotalc[i] += timec[i];
            }

            time = timec[0];
            if( ncu>0 )time += (timec[1]/dncu)*(ncu/ncf)*dncf;
            if( ncv>0 )time += (timec[2]/dncv)*(ncv/ncf)*dncf;
            
            previous_ncv=ncv;
            previous_ncu=ncu;
            previous_ncf=ncf;
         }

         break;
      case -1:
         // for reproducible debugging tests
         if(ncalls==0)srand( getRank(0) );
         ncalls++;

         time = 10.+(double)rand()/(double)RAND_MAX;
         break;
   }
   
   return time;   
}

PTIMING *profile(int field, int FLAG)
{
	double alpha;
   switch (FLAG)
   {
   case INIT:
      for (int i = 0; i < NTIMING; i++)
      {
         timing[i].status=OFF; 
         timing[i].elapsetime = timing[i].elapsetotal = 0.0;
         timing[i].elapseave = -2;
         timing[i].cputime = timing[i].cputotal = 0.0;
         timing[i].cpuave = 0.0;
         timing[i].ncalls = 0.0;
      }
      timing[TOTAL].name = "Total";
      timing[LOOP].name = "Loop";
      timing[MDSTEP].name = "   MDSTEP";
      timing[GROUPS].name = "   GROUPS";
      timing[KINETIC_TERMS].name = "     DDCKin";
      timing[LOOPIO].name = "   LoopIO";

      timing[ASSIGNMENT].name    = "Assignment";
      timing[ASSIGNMENT_B0].name = "   Assign_b0";
      timing[DOMAIN_SET].name    = "      domSet";
      timing[ASSIGNMENT_B2].name = "   Assign_b2";
      timing[ASSIGNMENT_B3].name = "   Assign_b3";
      timing[B3_1].name          = "      b3.1";
      timing[B3_2].name          = "      b3.2";
      timing[B3_3].name          = "      b3.3";
      timing[B3_4].name          = "      b3.4";
      timing[ASSIGNMENT_B5].name = "   Assign_b5";
      timing[COMMLIST].name    = "CommList";
      timing[COMMDOMS].name    = "   DomainSz";
      timing[COMMALLG].name    = "     CommAllg";
      timing[COMMCALC].name    = "   CommCalc";
      timing[FIND_REMOTE].name = "     FindRmte";
      timing[COMMCOMM].name    = "   CommComm";
      timing[FINDPAIR].name    = "   FindPair";
      timing[WAIT_FPAIRS].name = " WaitFPairs";
      timing[BOX].name         = "    GeomBox";
      timing[BOX_FAST].name    = "     GmFast";
      timing[BOX_DEFAULT].name = "     GmDflt";
      timing[UPDATE].name = "Update";
      timing[UPDATEWAIT].name = "  UpdateWait";
      timing[P_UPDATEFORCE].name = "UpdateForc";
      timing[BARRIER1].name = "Barrier1";

      timing[DDCENERGY].name = "DDCEnergy";
      timing[UPDATEALL].name = "  UpdateAll";
      timing[UPDATETABLES].name = "  UpdatTab";

      timing[PARITYCHECK].name = "ParityChk";
      timing[SAVESTATE].name = "SaveState";

      timing[P_FORCE].name      	     = "Force";
      timing[CUTOFF].name             = "Cutoff"; 
//      timing[MEMREAD].name            = "  memRead";
      timing[P_KER].name     	        = "  Kernel";
      timing[TFBIN].name              = "  TFBin"; 
      timing[TFMERGE].name            = "  TFMerge"; 
      timing[TFZBAR].name             = "  TFzBar"; 
      timing[TFSOLVE].name            = "  TFSolve"; 
      timing[TFPKAPPAZBAR].name       = "  TFKappazBar"; 
      timing[MGPT_MANYBODY].name      = "  MGPT_ManyBody";
      timing[CHARMM_T].name           = "   CHARMM_Tot"; 
      timing[CHARMM_COVALENT].name    = "    COVALENT"; 
      timing[CHARMM_CONNECTIVE].name  = "     CONNECT"; 
      timing[CHARMM_MALLOC].name      = "    MALLOC"; 
      timing[CHARMM_NONBOND].name     = "    NONBOND"; 
      timing[EWALD_T].name 		     = "   ewald";
      timing[RENDEZVOUS_T].name 	     = "    RonDayVu";
      timing[KSPACE_BEGIN].name 	     = "    kspaceB";
      timing[RHOMAKE].name 		     = "     rho_mak";
      timing[GRIDDEN].name 		     = "     gridDen";

      timing[GATHER_REDUCE].name      = "     GathRed";
      timing[GATHER_REDUCELOOP1].name = "      GRLOOP1";
      timing[GATHER_REDUCELOOP2].name = "      GRLOOP2";
      timing[EGR_BARRIER].name        = "     egr_bar";
      timing[PGATHER].name            = "     Pgather";
      timing[POISSON].name 		     = "     poisson";
      timing[FFT_FORWARD].name 	     = "      fft_f";
      timing[FFT_FORWARD_COMM].name   = "       fcomm";
      timing[FFT_BACKWARD].name	     = "      fft_b";
      timing[FFT_BACKWARD_COMM].name  = "       bcomm";
      timing[PSCATTER].name 	        = "     Pscatte";
      timing[EKB_BARRIER].name	     = "    ekb_barr";
      timing[RSPACE].name 		        = "    rspace";
      timing[BUILD_SCHEDULE].name     = "   schedule";
      timing[PAIREVAL].name     	     = "   pairEval";
      timing[KSPACE_END].name 	     = "   kspaceE";
      timing[UNPACK_GRID].name 	     = "    unpack";
      timing[UNPACK_WAIT].name        = "     wait";
      timing[SCATTER_POTENTIAL].name  = "    ScatPot";
      timing[SCATTER_POTLOOP1].name   = "      SPLOOP1";
      timing[SCATTER_POTLOOP2].name   = "      SPLOOP2";
      timing[SS_WAIT].name 		     = "     ssWait";
      timing[ELECTRIC].name 		     = "    electri";
      timing[OPT_TRANS].name 		     = "    opt_trans";
      timing[P_OPT_BOX].name 		     = "    opt_box";
      timing[BLOCK0].name             = "    block0";
      timing[EAM_PASS1].name 		     = "    eamPass1";
      timing[EAM_EMBED].name 		     = "    eamEmbed";
      timing[EAM_PASS2].name 		     = "    eamPass2";
      timing[VERLET].name     	     = "Verlet";
      timing[VBOX].name       	     = "  GeomBox";
      timing[PAIRLIST].name   	     = "  PairList";
      timing[UPDATEPAIR].name 	     = "UpdatePair";
      timing[EVAL_ETYPE].name         = "energyInfo";
      timing[BARRIER2].name           = "Barrier2";
      timing[BARRIER3].name           = "Barrier3";
      timing[LDBAL].name              = "LoadBalance";
      timing[QHULL].name              = "   QHull";


		/* Timers for recursive bisection code */ {
		  const int idx[] = {
			 RECBIS_PARKSTAT  ,RECBIS_KSTAT  ,RECBIS_BALANCE  ,RECBIS_EXCHANGE  ,RECBIS_EQUALIZE,
			 RECBIS_BALANCE2  ,RECBIS_A2ASPARSE,RECBIS_REDIST,RECBIS_REDIST2CORE,RECBIS_BACKCOMM
		  };
		  char *names[] = {
			 "RECBIS_PARKSTAT   ",
			 "RECBIS_KSTAT      ",
			 "RECBIS_BALANCE    ",
			 "RECBIS_EXCHANGE   ",
			 "RECBIS_EQUALIZE   ",
			 "RECBIS_BALANCE2   ",
			 "RECBIS_A2ASPARSE  ",
			 "RECBIS_REDIST     ",
			 "RECBIS_REDIST2CORE",
			 "RECBIS_BACKCOMM   "
		  };
		  for(unsigned i = 0; i<sizeof(idx)/sizeof(*idx); i++)
			 timing[idx[i]].name = names[i];
		}


      /* reset counters */
      profileN(-2);
      profileN2Search(-2);
      profileN2(-2);
      profileN2R(-2);
      profileMT(-2);
      profileNBox(-2);
      profileThrdImbalance(-2);
      break;
   case ZERO:
      for (int i = 0; i < NTIMING; i++)
		{
			if (i != TOTAL)
			{
				timing[i].elapsetime = timing[i].elapsetotal = 0.0;
				timing[i].elapseave = -2;
				timing[i].cputime = timing[i].cputotal = 0.0;
				timing[i].cpuave = 0.0;
				timing[i].ncalls = 0.0;
			}
		}
		break;
	case ACCUM:
		timing[field].status=ON; 
		timing[field].elapsetime = 0.0;
		timing[field].cputime = 0.0;
		timing[field].elapsestart = pelapsetime();
		timing[field].cpustart = pcputime();

		break;
	case START:
		{
			timing[field].status=ON; 
			timing[field].ncalls += 1.0;
			timing[field].elapsetime = 0.0;
			timing[field].cputime = 0.0;
			timing[field].elapsestart = pelapsetime();
			timing[field].cpustart = pcputime();
/*
			if (field == P_FORCE)
			{
				printf("FORCE START %d %f\n",P_FORCE,field,timing[field].elapsestart);
				fflush(stdout); 
			}
*/
		}
		break;
	case END:
		timing[field].status=OFF; 
		double et = pelapsetime() - timing[field].elapsestart;
		double ct = pcputime() - timing[field].cpustart;
		timing[field].elapsetime += et;
		timing[field].cputime += ct;
		timing[field].elapsetotal += et;
		timing[field].cputotal += ct;
/*
		if (field == P_FORCE)
		{
			printf("FORCE %d %d %f %f\n",P_FORCE,field,et,timing[field].elapsetotal);
			fflush(stdout); 
		}
*/

		break;
	case AVE:
      for (int i = 0; i < NTIMING; i++)
		{
			alpha = 0.5;
			if (timing[i].elapseave > -1)
			{
				timing[i].elapseave = (1.0 - alpha)*timing[i].elapseave + alpha*timing[i].elapsetime;
				timing[i].cpuave = (1.0 - alpha)*timing[i].cpuave + alpha*timing[i].cputime;
			}
			else
			{
				timing[i].elapseave = timing[i].elapsetime;
				timing[i].cpuave = timing[i].cputime;
			}
		}
		break;
	}
	fflush(stdout);
	return timing;
}

void profileWrite(PFILE*file, int id)
{
	double e, et, c, ct, ea, ca;
	double     ets, cts;
	double Ces, Ccs, Ceas, Ccas, Cets, Ccts;
	double Ees, Ecs, Eeas, Ecas, Eets, Ects;
	char *s;
   double naverage = profileN(-1);
   double n2Searchaverage = profileN2Search(-1);
   double n2average = profileN2(-1);
   double n2Raverage = profileN2R(-1);
   double mtaverage = profileMT(-1);
   double imbaverage= profileThrdImbalance(-1);
   int    numSchedule = profileMT(-3);
 
	Pprintf(file, "Timings for Processor %d", id);
	if (hi_hasTorus() != 0)
	{
      int nDim = hi_nTorusDim();
      int coord[nDim];
      hi_torusCoords(coord);
      int coreId = hi_coreId();
	   Pprintf(file, "  Torus Coords (coreId last) = %d", coord[0]);
      for (int ii=1; ii<nDim; ++ii)
         Pprintf(file, " %d", coord[ii]);
      Pprintf(file, " %d\n", coreId);
	}
// double es = 0.0; 
	ets = cts = 0.0;
	Ces = Ccs = Ceas = Ccas = Cets = Ccts = 0.0;
	Ees = Ecs = Eeas = Ecas = Eets = Ects = 0.0;
   const double eloop = timing[LOOP].elapsetotal;
   for (int i = 0; i < NTIMING; i++)
	{
		e = timing[i].elapsetime;
		c = timing[i].cputime;
		et = timing[i].elapsetotal;
		ct = timing[i].cputotal;
		ea = timing[i].elapseave;
		ca = timing[i].cpuave;
		s= timing[i].name;
		if ((i == ASSIGNMENT) || (i == COMMLIST) || (i == UPDATE) || (i == P_UPDATEFORCE) || (i == BARRIER1) )  //  Communication Total 
		{
			Cets += et;
			Ccts += ct;
			Ces += e;
			Ccs += c;
			Ceas += ea;
			Ccas += ca;
		}
		if ((i == P_FORCE) || i == CUTOFF || (i == VERLET) || i== UPDATEPAIR || i==EVAL_ETYPE || i == BARRIER2 || i == BARRIER3)     // Energy Calc Total 
		{
			Eets += et;
			Ects += ct;
			Ees += e;
			Ecs += c;
			Eeas += ea;
			Ecas += ca;
		}
		ets = Cets + Eets;
		cts = Ccts + Ects;
		//es = Ces + Ees;
		//eas = Ceas + Eeas;
		//cas = Ccas + Ecas;
   }
	Pprintf(file, "\n");
	Pprintf(file, "------------------------------------------------------------------------\n");
	Pprintf(file, "                        |             Elaspe Timings                 \n");
	Pprintf(file, "                 #Calls |  Current  Average  Aging  Total     %%Loop   \n");
	Pprintf(file, "------------------------------------------------------------------------\n");
   /*
		if (i == DDCENERGY) 	Pprintf(file, "------------------------------------------------------------------------\n");
		if (i == PARITYCHECK) 	Pprintf(file, "------------------------------------------------------------------------\n");
		if (i == ASSIGNMENT) 	Pprintf(file, "------------------------------------------------------------------------\n");
		if (i == LDBAL) 	Pprintf(file, "------------------------------------------------------------------------\n");
      */
	Pprintf(file, "------------------------------------------------------------------------\n");
   for (int i = TOTAL;i<ASSIGNMENT; i++)
   {
      e = timing[i].elapsetime;
      c = timing[i].cputime;
      et = timing[i].elapsetotal;
      ct = timing[i].cputotal;
      ea = timing[i].elapseave;
      ca = timing[i].cpuave;
      double nC_i = timing[i].ncalls;
      s= timing[i].name;
      if (timing[i].ncalls > 0)
      {
            if (i == DDCENERGY || i ==LOOP) 	Pprintf(file, "------------------------------------------------------------------------\n");
            Pprintf(file, 			"%-13s: %7.0f  | %7.2f %7.2f %7.2f %7.2f %10.2f \n", s, nC_i, e, et/nC_i, ea, et, et/eloop*100.0);
      }
   }
   Pprintf(file, "------------------------------------------------------------------------\n");
   Pprintf(file, "------------------------------------------------------------------------\n");

   for (int i = ASSIGNMENT; i<= BARRIER1; i++)
   {
      e = timing[i].elapsetime;
      c = timing[i].cputime;
      et = timing[i].elapsetotal;
      ct = timing[i].cputotal;
      ea = timing[i].elapseave;
      ca = timing[i].cpuave;
      double nC_i = timing[i].ncalls;
      s= timing[i].name;
      if (timing[i].ncalls > 0)
         Pprintf(file, 			"%-13s: %7.0f  | %7.2f %7.2f %7.2f %7.2f %10.2f \n", s, nC_i, e, et/nC_i, ea, et, et/eloop*100.0);
   }
   double nC = timing[MDSTEP].ncalls;
   Pprintf(file, "------------------------------------------------------------------------\n");
   Pprintf(file, 				"%-13s: %7.0f  | %7.2f %7.2f %7.2f %7.2f %10.2f \n",
         "Comm Total", nC, Ces, Cets/nC, Ceas, Cets, Cets/eloop*100.0);

   Pprintf(file, "------------------------------------------------------------------------\n");
   Pprintf(file, "------------------------------------------------------------------------\n");

   for (int i = P_FORCE; i <= BARRIER3; i++)
   {
      e = timing[i].elapsetime;
      c = timing[i].cputime;
      et = timing[i].elapsetotal;
      ct = timing[i].cputotal;
      ea = timing[i].elapseave;
      ca = timing[i].cpuave;
      double nC_i = timing[i].ncalls;
      s= timing[i].name;
      if (timing[i].ncalls > 0)
         Pprintf(file, 			"%-13s: %7.0f  | %7.2f %7.2f %7.2f %7.2f %10.2f \n", s, nC_i, e, et/nC_i, ea, et, et/eloop*100.0);
   }
   Pprintf(file, "------------------------------------------------------------------------\n");
   Pprintf(file, 				"%-13s: %7.0f  | %7.2f %7.2f %7.2f %7.2f %10.2f \n",
         "Energy Total", nC, Ees, Eets/nC, Eeas, Eets, Eets/eloop*100.0);
   Pprintf(file, "------------------------------------------------------------------------\n");
   Pprintf(file, "------------------------------------------------------------------------\n");

   Pprintf(file, "Average # Particles           = %f\n", naverage);
   Pprintf(file, "Average # pairSearched        = %f\n", n2Searchaverage);
   Pprintf(file, "Average # pairsInList         = %f\n", n2average);
   Pprintf(file, "Average # pairsCalculated     = %f\n", n2Raverage);
   Pprintf(file, "Average # Boxes               = %f\n", profileNBox(-1));
   Pprintf(file, "Thread schedule method        = %s\n", profileThreadScheduleMethodName(NULL));
   Pprintf(file, "Average thread imbalance      = %f\n", imbaverage);
   Pprintf(file, "Average Sum of Thread Memory  = %f\n", mtaverage);
   Pprintf(file, "# Thread Schedule             = %d\n", numSchedule);
   Pprintf(file, "Force (Elapse Time)/(particles *steps)  nsec = %.3e\n", 1e9*timing[P_FORCE].elapsetotal/(naverage*timing[P_FORCE].ncalls));
   Pprintf(file, "Total (Elapse Time)/(particles*steps) nsec = %.3e\n", 1e9*ets/(naverage*timing[P_FORCE].ncalls));
   Pprintf(file, "MDStep  (Elapse Time)/(particles *steps)  nsec = %.3e\n", 1e9*timing[MDSTEP].elapsetotal/(naverage*timing[P_FORCE].ncalls));
   Pprintf(file, "Loop  (Elapse Time)/(particles *steps)  nsec = %.3e\n", 1e9*timing[LOOP].elapsetotal/(naverage*timing[P_FORCE].ncalls));
   {
      extern int nresizings,buf1cap,buf2cap;
      Pprintf(file, "@ Number of resizings       = %f\n", (double) nresizings);
      Pprintf(file, "@ buf1Capacity              = %f\n", (double) buf1cap);
      Pprintf(file, "@ buf2Capacity              = %f\n", (double) buf2cap);
      /* Dump data from rectimer (e.g. recursive bisection functions) */ 
      {
         if(0) {
            char *str = rectimer_printreport_string();
            //if(getRank(0) == 0) rectimer_printreport(stdout);
            Pprintf(file, "@ %s @ -----\n\n",str);
            free(str);
         } else {
			  rectimer_printreport_file(Pprintf,file);
         }
      }
   }

   Pprintf(file, "\n");
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
