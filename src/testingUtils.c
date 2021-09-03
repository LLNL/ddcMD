#include "commandLineOptions.h"
#include "testingUtils.h"


void objectSetup(COMMAND_LINE_OPTIONS  opt, MPI_Comm comm);
void mpiStartUp(int argc, char *argv[]);
void   parityErrorReport();

void version_init(int argc,char *argv[]);
void simulate(COMMAND_LINE_OPTIONS opt, MPI_Comm SimulateComm);
void commons_init(void);

extern MPI_Comm COMM_LOCAL;

static ROUTINE *routine;
//static char* testFile = "";


//
char *exename; // just a hack to get code to compile.  Shiv will need to make the real fix. 
SIMULATE* setupSimulate( int i) {
    
    int argc = 1;
    char **argv = (char**) malloc(2*sizeof(char*));
    argv[0]=exename;;
    argv[1]="\0";
   COMMAND_LINE_OPTIONS opt = parseCommandLine(argc, argv);
   //if (i==0){
     mpiStartUp(argc,argv);
   //}
   checkLimits();
   commons_init();
   prime_init(30000, getRank(-1), getSize(-1));
   units_internal(a0_MKS,Rinfhc_MKS*1e-30/(a0_MKS*a0_MKS),1e-15,e_MKS/1e-15,Rinfhc_eV/kB_eV,1.0,1.0);
   units_external(1e-10,u_MKS,1e-15,e_MKS/1e-15,1.0,1.0,1.0);

   setParityHandler();
   version_init(argc,argv);

        profile(NOFIELD, INIT);
   HPM_Init();

   profile(TOTAL, START);

   HPM_Start("Total");
   WAIT(-1);
   objectSetup(opt, MPI_COMM_WORLD);
   MASTER master=masterFactory(opt);
   //
   routine  = routineManager_init(NULL,"routineManager",master.fcn,master.parms);
   SIMULATEMASTERPARMS *smParms = (SIMULATEMASTERPARMS *)master.parms;
   char * name = (char*) malloc(16*sizeof(char));
   strcpy(name, (smParms->common.simulateName));
   printf("name: %s\n", name);
   if (i>1)
   {  
      strcpy(name,smParms->common.simulateName);
      printf("%s%i",name, i);
      smParms->common.simulateName = (char*) &name;
   }
   SIMULATE* simulate = simulate_init(NULL,smParms->common.simulateName, MPI_COMM_WORLD);
   adjustBox(simulate);
   return simulate;
}


/*
SIMULATE* setupSimulate2() {
    int argc = 1;
    char **argv = (char**) malloc(2*sizeof(char*));
    argv[0]=exename;;
    argv[1]="\0";
   COMMAND_LINE_OPTIONS opt = parseCommandLine(argc, argv);
   checkLimits();
   commons_init();
   prime_init(30000, getRank(-1), getSize(-1));
   units_internal(a0_MKS,Rinfhc_MKS*1e-30/(a0_MKS*a0_MKS),1e-15,e_MKS/1e-15,Rinfhc_eV/kB_eV,1.0,1.0);
   units_external(1e-10,u_MKS,1e-15,e_MKS/1e-15,1.0,1.0,1.0);

   setParityHandler();
   version_init(argc,argv);

        profile(NOFIELD, INIT);
   HPM_Init();

   profile(TOTAL, START);

   HPM_Start("Total");
   WAIT(-1);
   objectSetup(opt, MPI_COMM_WORLD);
   SIMULATE* simulate = simulate_init(NULL,"simulate2", MPI_COMM_WORLD);
   adjustBox(simulate);
   return simulate;
}
*/






void setupFirstEnergyCall(SIMULATE * simulate){
   SYSTEM * sys = simulate->system;
   COLLECTION *collection = sys->collection;
   STATE *state = collection->state; 
   ddc_put(simulate->ddc, PAIRMODE, LABEL);
   ddc_put(simulate->ddc, CORNER, &sys->box->corner);
   ddc_put(simulate->ddc, DDCNLOCAL, sys->nlocal);
   simulate->ddc->update = 3;
    
   resize(collection->size, 2, collection->state); // ensure that any quantities registered with particleRegisterinfo are allocated before we try to use them in ddcenergy.
   for (int kk=0; kk<sys->ngroup; kk++) 
   {   
      sys->group[kk]->start(sys->group[kk],START_GROUP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
      sys->group[kk]->Update1(sys->group[kk],START_GROUP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
   }   

   timestamp("Starting First Call to ddcenergy");

}
//
void setupEnergyCall(SIMULATE * simulate){
	int e_eval_flag = 1; 
	SYSTEM * sys = simulate->system;
	DDC * ddc = simulate->ddc;
        STATE *state=sys->collection->state;
        profile(DDCENERGY, START);
        if (sys->potential[0] ->itype == ZEROPOTENTIAL && sys->npotential==0 )   // specical case of no interaction potential. Only kinetic terms needed. 
        {
                zeroAll(sys);
                if (e_eval_flag)
                {
                        kinetic_terms(sys, 1);
                        eval_energyInfo(sys);
                }
                return ;
        }
        //ETYPE *energyInfo = &sys->energyInfo;

        cutoffs(ddc,sys);
        auxNeighbor_begin(sys->nion);
        int update=0;
        if (e_eval_flag == -1)
        {
                e_eval_flag =0;
                update=-1;
        }
        {   
                if ( TEST0(sys->loop, ddc->loadBalance->rate) )
                {   
                        ddc->loadBalance->balanceFunction(ddc->loadBalance, sys, ddc);
                        update=1;
                }   
        }   
    
   profile(UPDATEALL, START);
//        switch (sys->potential[0]->itype)
//        {   
//                case EAM_OPT:
//                case EAM_ONEPASS:
//                case PAIRENERGY:
//                       ddcUpdateAll_pair(ddc, sys, energyInfo,update);  
//                        break; 
//                case EAM1PASS:
//                case EAM2PASS:
//                        ddcUpdateAll_(ddc, sys, energyInfo, update);  
//                        break; 
//                case PLASMA:
//                        if (sys->potential[0]->call_fsumX) ddcUpdateAll(ddc, sys, energyInfo, update);  
////                      else ddcUpdateAll_pair(ddc, sys, energyInfo, update);  
//                        else ddcUpdateAll_(ddc, sys, energyInfo, update);  
//                        break; 
//                case EWALD:
//                        if (sys->potential[0]->call_fsumX) ddcUpdateAll(ddc, sys, energyInfo, update);  
//                        break; 
//              default:
//                        ddcUpdateAll(ddc, sys, energyInfo, update);  
//        }   
        profile(UPDATEALL, END);
        zeroAll(sys);
        profile(BARRIER1, START);
        WAIT(0);
        profile(BARRIER1, END);
        profile(P_FORCE, START);
   system_pCalculate(sys); 
        if (state->q != NULL) for (unsigned i=0;i<sys->nion;i++) state->q[i] = ((ATOMTYPE_PARMS*)(state->species[i])->parm)->charge; 

	//energy call	
        //for (int i = 0; i < sys->npotential; i++) sys->potential[i]->eval_potential((void *)sys, (void *)sys->potential[i]->parms, &(sys->energyInfo)); 

	//stuff that goes after energy call. this is now in finishEnergyCall()
/*
        if (sys->potential[0]->call_fsumX) fsumX(sys, &(sys->energyInfo.virial));
        profile(P_FORCE, END);
        profile(BARRIER2, START);
    //    idleThreadsCheck(sys, COMM_LOCAL);  
     //   if (parityCheckingBarrier(COMM_LOCAL) != 0) return -1; 
        profile(BARRIER2, END);
        ddcUpdateForce(ddc,sys->pCalculate);  
        //ddcUpdateAccum(ddc); 
        if (sys->energyInfo.pdiaV != 0.0) for (unsigned i=0;i<sys->nlocal;i++)  state->virial[i]+=sys->energyInfo.pdiaV/sys->nglobal; 

        profile(BARRIER3, START);
        WAIT(0);
        profile(BARRIER3, END);
        if (e_eval_flag) 
        {   
                kinetic_terms(sys, 1); 
                eval_energyInfo(sys);
        }   
        profile(DDCENERGY, END);
//        return 0;
*/
}
//
void finishEnergyCall(SIMULATE * simulate){
	int e_eval_flag =1 ;
	SYSTEM * sys = simulate->system; 
	DDC * ddc = simulate->ddc; 
	STATE * state = sys->collection->state;

        if (sys->potential[0]->call_fsumX) fsumX(sys, &(sys->energyInfo.virial));
        profile(P_FORCE, END);
        profile(BARRIER2, START);
        profile(BARRIER2, START);
    //    idleThreadsCheck(sys, COMM_LOCAL);  
    //         //   if (parityCheckingBarrier(COMM_LOCAL) != 0) return -1; 
        profile(BARRIER2, END);
        ddcUpdateForce(ddc,sys->pCalculate);
        if (sys->energyInfo.pdiaV != 0.0) for (unsigned i=0;i<sys->nlocal;i++)  state->virial[i]+=sys->energyInfo.pdiaV/sys->nglobal;

        profile(BARRIER3, START);
        WAIT(0);
        profile(BARRIER3, END);
        if (e_eval_flag)
        {
                kinetic_terms(sys, 1);
                eval_energyInfo(sys);
        }
        profile(DDCENERGY, END);
}

/*
void setupEnergyCall(SIMULATE * simulate)
{
	int e_eval_flag=1;
	DDC * ddc = simulate->ddc;
	SYSTEM *sys = simulate->system;
	// start ddcenergy
        STATE *state=sys->collection->state;
        profile(DDCENERGY, START);
        if (sys->potential[0] ->itype == ZEROPOTENTIAL && sys->npotential==0 )   // specical case of no interaction potential. Only kinetic terms needed. 
        {
                zeroAll(sys);
                if (e_eval_flag)
                {
                        kinetic_terms(sys, 1);
                        eval_energyInfo(sys);
                }
                return 0;
        }
        ETYPE *energyInfo = &sys->energyInfo;

        cutoffs(ddc,sys);
        auxNeighbor_begin(sys->nion);
        int update=1;
        if (e_eval_flag == -1)
        {
                e_eval_flag =0;
                update=-1;
        }

        //some comment goes here
        {
                SIMULATE* simulate = simulate_getSimulate(NULL);
                if ( TEST0(sys->loop, ddc->loadBalance->rate) )
                {
                        ddc->loadBalance->balanceFunction(ddc->loadBalance, sys, ddc);
                        update=1;
                }
        }

   profile(UPDATEALL, START);
        switch (sys->potential[0]->itype)
        {
                case EAM_OPT:
                case EAM_ONEPASS:
                case PAIRENERGY:
                        ddcUpdateAll_pair(ddc, sys, energyInfo,update);
                        break;
                case EAM1PASS:
                case EAM2PASS:
                        ddcUpdateAll_(ddc, sys, energyInfo, update);
                        break;
                case PLASMA:
                        if (sys->potential[0]->call_fsumX) ddcUpdateAll(ddc, sys, energyInfo, update);
                        else ddcUpdateAll_pair(ddc, sys, energyInfo, update);
                        break;
                case EWALD:
                        if (sys->potential[0]->call_fsumX) ddcUpdateAll(ddc, sys, energyInfo, update);
                        break;
                default:
                        ddcUpdateAll(ddc, sys, energyInfo, update);
        }
        profile(UPDATEALL, END);
        zeroAll(sys);

        profile(BARRIER1, START);
        WAIT(0);
        profile(BARRIER1, END);
        profile(P_FORCE, START);
   system_pCalculate(sys);
        if (state->q != NULL) for (unsigned i=0;i<sys->nion;i++) state->q[i] = ((ATOMTYPE_PARMS*)(state->species[i])->parm)->charge;



        //call ewald cpu



////call ewald gpu



printf("done\n");
}
*/
