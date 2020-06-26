#include <string.h>

#include "simulate.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "bioCharmmCovalent.h"
#include "bioGid.h"
#include "ddcRule.h"
#include "ddc.h"
#include "ddcMalloc.h"
#include "utilities.h"
char *strcat_(char *a, char *b);

void ddcRuleCharmm_init(DDCRULE *ddcRule)
{
    DDC* ddc = (DDC*) ddcRule->parent;
    SIMULATE* simulate = (SIMULATE*) (ddc->parent);
    SYSTEM* system = simulate->system;
    int foundCHARMM_PARMS = 0;
    for (int i = 0; i < system->npotential; ++i)
    {
        //if (strcmp(system->potential[i]->type, "CHARMM") == 0){
        if (system->potential[i]->itype == CHARMM)
        {
            CHARMMPOT_PARMS* charmmpot_parms = (CHARMMPOT_PARMS *) (system->potential[i]->parms);
            CHARMM_PARMS *charmmParms = charmmpot_parms->charmmParms;
            ddcRule->parms = charmmParms;
            foundCHARMM_PARMS = 1;
            break;
        }
    }
    if (foundCHARMM_PARMS == 0)
    {
        error_action("ddcRule Error: Cannot find CHARMM_PARMS in the system potentials", ERROR_IN("ddcRule", ABORT));
    }
}

void ddcRuleCharmm(DDC* ddc)
{

    //MPI_Comm comm = MPI_COMM_WORLD;
    //int myrank;
    //MPI_Comm_rank(comm, &myrank);

    DDCRULE* ddcRule = ddc->ddcRule;
    CHARMM_PARMS *charmmParms = (CHARMM_PARMS *) ddcRule->parms;

    PARTICLE* particles = ddc->particles;
    int nlocal = ddc->number_local;
    int nremote = ddc->number_remote;
    int nion = nlocal + nremote;

    //for (int i = 0; i < nion; i++) {
    //    printf("rank=%d Particle i=%d gid=%llu doID=%d\n", myrank, i, particles[i].global_index, particles[i].domain_id);
    //}

    GID_ORDER* gidList = (GID_ORDER*) ddcMalloc(nion * sizeof (GID_ORDER));
    for (int i = 0; i < nion; i++)
    {
        gidList[i].id = i;
        gidList[i].gid = particles[i].global_index;
        //if(i<nlocal){
        //    printf("rank=%d LOCAL id=%d gid=%llu\n", myrank, gidList[i].id, gidList[i].gid);
        //}else{
        //    printf("rank=%d REMOTE id=%d gid=%llu\n", myrank, gidList[i].id, gidList[i].gid);
        //}
    }

    qsort(gidList, nion, sizeof (GID_ORDER), compare);

    //for (int i = 0; i < nion; i++) {
    //  printf("qsort rank=%d qsort id=%d gid=%llu\n", myrank, gidList[i].id, gidList[i].gid);
    //}

    int nalloc = 100;
    RESI_CONN** resiConnList = ddcMalloc(sizeof (RESI_CONN*) * nalloc);
    int rlist = 0;

    gid_type old_molres = ~0; // Make sure it is different at beginning

    //Find list of residues
    int totNumPad = 0;
    RESI_CONN* resiConn = NULL;
    //char specieName[11];

    //char* oldName=specieName;
    for (int i = 0; i < nion; i++)
    {
        gid_type gid = gidList[i].gid;
        gid_type molres = gid&molresMask;
        if (molres != old_molres)
        {
            if (rlist >= nalloc)
            {
                nalloc = nalloc * 2 + 100;
                resiConnList = (RESI_CONN**) ddcRealloc(resiConnList, sizeof (RESI_CONN*) * nalloc);
            }
            int index = particles[gidList[i].id].type;
            SPECIES* species = species_by_index(NULL, index & 0xffff);

            //if(strcmp(species->name, oldName)!=0){
            //Don't want to repeatedly search for the same name
            resiConn = findResiConnNew(charmmParms, species->name);
            //strcpy(oldName, species->name);
            //}
            //printf("species->name=%s, index=%d\n", species->name, index);
            resiConnList[rlist] = resiConn;
            totNumPad += resiConn->atomListSize;
            old_molres = molres;
            rlist++;
        }

    }

    //for (int i = 0; i < rlist; i++) {
    //    printf("Rank=%d i=%d name=%s numAtom=%d numHeavy=%d numHydrogen=%d\n", 
    //            myrank, i, resiConnList[i]->resName, resiConnList[i]->atomListSize,    //resClean
    //            resiConnList[i]->heavyAtmsSize, resiConnList[i]->heavyAtms[0]->numHydrogen);
    //}

    //determine the residue range
    RESRANGE resRange[rlist];
    unsigned start = 0;
    for (int i = 0; i < rlist; i++)
    {
        resRange[i].start = start;
        start += resiConnList[i]->atomListSize;
        resRange[i].end = start;
    }

    //printf("Padding\n");
    //padding
    //PARTICLE* particlepad = (PARTICLE*) ddcMalloc(totNumPad * sizeof (PARTICLE));
    int *exist = (int *) ddcMalloc(totNumPad * sizeof (int));
    int oldIp = -1;
    int ip = 0;
    int ippad = 0;
    SPECIES* species = NULL;
    //SPECIES* species = species_by_index(NULL, particles[gidList[ip].id].type & 0xffff);
    for (int i = 0; i < rlist; i++)
    {
        RESI_CONN* resiConn = resiConnList[i];
        char * name = resiConn->resName;
        char *terminus = "x";
        if (resiConn->nTer)
        {
            terminus = "n";
        }
        if (resiConn->cTer)
        {
            terminus = "c";
        }
        for (int j = 0; j < resiConn->atomListSize; j++)
        {
            char atm[11];
            char *atmName = atm;
            strcpy(atmName, name);
            strcat_(atmName, terminus);
            strcat_(atmName, resiConn->atomList[j]->atmName);

            //printf("rank=%d atmName=%s size=%d species->name=%s sp_size=%d\n", 
            //        myrank, atmName, strlen(atmName), species->name, strlen(species->name));
            //int atmSize=strlen(atmName);
            //int spSize=strlen(species->name);
            if (ip != oldIp)
            {
                species = species_by_index(NULL, particles[gidList[ip].id].type & 0xffff);
                oldIp = ip;
            }

            if (strcmp(atmName, species->name) == 0)
            {
                //particlepad[ippad] = particles[gidList[ip].id];
                //particlepad[ippad].domain_id = particles[gidList[ip].id].domain_id;

                //printf("rank=%d ip=%d gidList[ip].id=%d atmName=%s species->name=%s\n", 
                //        myrank, ip, gidList[ip].id, atmName, species->name);

                exist[ippad] = gidList[ip].id;

                ip++;
                //species = species_by_index(NULL, particles[gidList[ip].id].type & 0xffff);

            }
            else
            {
                exist[ippad] = -1;
            }
            ippad++;
        }
    }

    //   for (int i = 0; i < totNumPad; i++) {
    //       printf("rank=%d exist[%d]=%d\n", myrank, i, exist[i]);
    //   }

    //MPI_Barrier(comm);
    for (int i = 0; i < rlist; i++)
    {
        RESI_CONN* resiConn = resiConnList[i];

        if (resNameDiffer(resiConn->resName, "TIP3") == 0)
        {
            // Treat TIP3 water specially
            int j = 0;
            HEAVYATM* heavyAtm = resiConn->heavyAtms[j];
            int indexHeavyAtm = resRange[i].start + heavyAtm->atomID;
            if (exist[indexHeavyAtm] >= 0)
            {
                //printf("TIP3 rank=%d i=%d, resRange[i].start=%d, j=%d, atomID=%d indexHeavyAtm=%d\n",
                //        myrank, i, resRange[i].start, j, heavyAtm->atomID, indexHeavyAtm);
                for (int k = 0; k < heavyAtm->numHydrogen; k++)
                {
                    int indexH = resRange[i].start + heavyAtm->hydrogen[k];
                    if (exist[indexH] >= 0)
                    {
                        //flip the hydrogen domain id to be the same as its heavy atom
                        if (particles[exist[indexHeavyAtm]].domain_id != particles[exist[indexH]].domain_id)
                        {

                            //printf("rank=%d Hy_exist[%d]=%d, hy_doID=%d,H_exist[%d]=%d, h_doID=%d \n",
                            //        myrank, indexHeavyAtm, exist[indexHeavyAtm], particles[exist[indexHeavyAtm]].domain_id,
                            //        indexH, exist[indexH], particles[exist[indexH]].domain_id);
                            particles[exist[indexH]].domain_id = particles[exist[indexHeavyAtm]].domain_id;
                        }
                    }
                }
            }
        }
        else
        {
            // For the residue other than TIP3
            for (int j = 0; j < resiConn->heavyAtmsSize; j++)
            {
                HEAVYATM* heavyAtm = resiConn->heavyAtms[j];
                int indexHeavyAtm = resRange[i].start + heavyAtm->atomID;
                //printf("rank=%d i=%d, resRange[i].start=%d, j=%d, atomID=%d indexHeavyAtm=%d\n",
                //myrank, i, resRange[i].start, j, heavyAtm->atomID, indexHeavyAtm);
                if (exist[indexHeavyAtm] >= 0)
                {
                    // Heavy atom is in the local
                    for (int k = 0; k < heavyAtm->numHydrogen; k++)
                    {
                        int indexH = resRange[i].start + heavyAtm->hydrogen[k];
                        if (exist[indexH] >= 0)
                        {
                            //flip the hydrogen domain id to be the same as its heavy atom
                            if (particles[exist[indexHeavyAtm]].domain_id != particles[exist[indexH]].domain_id)
                            {

                                //printf("rank=%d Hy_exist[%d]=%d, hy_doID=%d,H_exist[%d]=%d, h_doID=%d \n",
                                //        myrank, indexHeavyAtm, exist[indexHeavyAtm], particles[exist[indexHeavyAtm]].domain_id,
                                //        indexH, exist[indexH], particles[exist[indexH]].domain_id);
                                particles[exist[indexH]].domain_id = particles[exist[indexHeavyAtm]].domain_id;
                            }
                        }
                    }

                }
            }
        }
    }

    //MPI_Barrier(comm);
    //for (int i = 0; i < nion; i++) {
    //    printf("END rank=%d Particle i=%d gid=%llu doID=%d\n", myrank, i, particles[i].global_index, particles[i].domain_id);
    //}

    //printf("END\n");
    ddcFree(exist);
    ddcFree(resiConnList);
    ddcFree(gidList);

}
