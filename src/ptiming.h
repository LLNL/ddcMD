#ifndef PTIMING_H
#define PTIMING_H

#ifdef __cplusplus
extern "C" {
#endif

struct pfile_st; // avoid including pio.h

enum PTIMING_ENUMS
  { NOFIELD, TOTAL, LOOP, MDSTEP,GROUPS, 
  KINETIC_TERMS, LOOPIO, DDCENERGY, UPDATEALL,UPDATETABLES,
  PARITYCHECK, SAVESTATE,
  LDBAL, QHULL,
  ASSIGNMENT, ASSIGNMENT_B0, DOMAIN_SET,
  ASSIGNMENT_B2, ASSIGNMENT_B3, B3_1, B3_2, B3_3, B3_4,
    ASSIGNMENT_B5, COMMLIST, COMMDOMS, COMMALLG, COMMCALC, FIND_REMOTE,
    COMMCOMM, FINDPAIR, WAIT_FPAIRS, BOX, BOX_FAST, BOX_DEFAULT, UPDATE, UPDATEWAIT, P_UPDATEFORCE,
  BARRIER1,
  P_FORCE, P_KER, CUTOFF, TFBIN,TFMERGE, TFZBAR,TFSOLVE, TFPKAPPAZBAR, MGPT_MANYBODY, 
  CHARMM_T, CHARMM_COVALENT,CHARMM_CONNECTIVE,  CHARMM_MALLOC, CHARMM_NONBOND,
  EWALD_T, RSPACE, RENDEZVOUS_T, KSPACE_BEGIN, RHOMAKE, GRIDDEN,
  GATHER_REDUCE, GATHER_REDUCELOOP1, GATHER_REDUCELOOP2, EGR_BARRIER,
  PGATHER, POISSON, FFT_FORWARD, FFT_FORWARD_COMM, FFT_BACKWARD, FFT_BACKWARD_COMM,
  PSCATTER, EKB_BARRIER, 
  BUILD_SCHEDULE,
  PAIREVAL, 
  KSPACE_END, 
  UNPACK_GRID, UNPACK_WAIT, SCATTER_POTENTIAL, SCATTER_POTLOOP1, SCATTER_POTLOOP2,SS_WAIT, ELECTRIC,
    OPT_TRANS, P_OPT_BOX, BLOCK0, EAM_PASS1, EAM_EMBED, EAM_PASS2,
  VERLET, VBOX, PAIRLIST, UPDATEPAIR,EVAL_ETYPE,BARRIER2,BARRIER3,
	 RECBIS_PARKSTAT,RECBIS_KSTAT,RECBIS_BALANCE,RECBIS_EXCHANGE,RECBIS_EQUALIZE,
	 RECBIS_BALANCE2,RECBIS_A2ASPARSE,RECBIS_REDIST,RECBIS_REDIST2CORE,RECBIS_BACKCOMM,
    GPU_BIN, GPUN_LIST, GPU_NO_LIST_PAIR, GPU_MALLOC, GPU_SEND, GPU_RECV,
  NTIMING};

enum PTIMING_ACTIONS { N, INIT, ZERO, ACCUM, START, END, AVE, WRITE, ON, OFF}; 

typedef struct ptiming_st
{
   double elapsestart;
   double cpustart;
   double elapsetime;
   double cputime;
   double elapsetotal;
   double cputotal;
   double elapseave;
   double cpuave;
   double ncalls;
   int status;
   char *name;
} PTIMING;

PTIMING *profile(int field, int FLAG);
PTIMING *profile_get_timing(void);
double profileN(int);
double profileN2Search(int n2);
double profileN2(int n2);
double profileN2R(int n2);
double profileNBox(int nbox);
double profileMT(int m);
char* profileThreadScheduleMethodName(char* st);
double profileThrdImbalance(double m);
void profileWrite(struct pfile_st* file, int myid);
double profileGetLoopTime(void); 
void profileSetLoopTime(double); 
void dumpprofile(char *filename); 
void slave_dumpprofile(char *filename); 
double profileGetLocalLoadElapsetime(int flag);

#ifdef __cplusplus
}
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
