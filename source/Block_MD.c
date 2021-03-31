//
//  Block_MD.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "Block_MD.h"

void Block_MD_St(void){

    int t, i, t0 = 0;
    double DT2 = (DT*DT), DToverMi, DTover2, DTover4, overMi, Mi, alpha;
    double cyclotronFreq, PARTMOM_Tilday, PARTPOS_Tilday;
    double sumLorentzForces, sumIntermolecularForces;
    struct point CF_t, SF_t, FLorentz;

    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};

    char therm;

    clock_t t_start, t_end;

    if (PRECONFIG_FLAG == 1) {

        FirstStep_St();
        t0 = 1;
    }
    for (t=t0; t<NTIMESTEPS; t++) {

        t_start = clock();

        alpha = 1.;
        therm = ' ';

        if (DEBUG_FLAG && _D_TOT_FORCES) CFtot.x = CFtot.y = CFtot.z = SFtot.x = SFtot.y = SFtot.z = Ftot.x = Ftot.y = Ftot.z = 0;

        if (THERMOSTAT == 'T' && fabs(ITEMP - TEMP) > (_TEMP_TOL*TEMP)){

            therm = '*';
            alpha = Thermostat(PARTPOS_TM1, PARTPOS_T, TEMP);

            if (DEBUG_FLAG && _D_THERMOSTAT) {

                printf("|∆T| = %.4e\t|∆T|_M = %.4e\talpha = %.4e\n", fabs(TEMP - ITEMP), _TEMP_TOL*TEMP, alpha);
            }
        }

        for (i=0; i<NPART; i++) {

            Mi = M[INDX[i]];
            overMi = 1./Mi;
            cyclotronFreq = Q[INDX[i]]*B0*overMi*0.5;
            DToverMi = DT*overMi;
            DTover2 = DT*0.5;
            DTover4 = DT*0.25;


            if (POT == 'J') {

                CF_t = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
                SF_t = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);

            } else if (POT == 'C') {

                CF_t = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
                SF_t = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
            }

            if (DEBUG_FLAG && _D_TOT_FORCES) {

                CFtot.x += (CF_t.x);
                CFtot.y += (CF_t.y);
                CFtot.z += (CF_t.z);

                SFtot.x += (SF_t.x);
                SFtot.y += (SF_t.y);
                SFtot.z += (SF_t.z);

                Ftot.x += (CF_t.x  + SF_t.x);
                Ftot.y += (CF_t.y  + SF_t.y);
                Ftot.z += (CF_t.z  + SF_t.z);
            }

            // collect forces
            sumIntermolecularForces += sqrt((CF_t.x + SF_t.x)*(CF_t.x + SF_t.x) + (CF_t.y + SF_t.y)*(CF_t.y + SF_t.y) + (CF_t.z + SF_t.z)*(CF_t.z + SF_t.z));
            FLorentz.x = Q[INDX[i]]*B0*PARTVEL[i].y;
            FLorentz.y = - Q[INDX[i]]*B0*PARTVEL[i].x;
            sumLorentzForces += sqrt(FLorentz.x*FLorentz.x + FLorentz.y*FLorentz.y);


            // SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + (PARTPOS_T[i].x - PARTPOS_TM1[i].x)*alpha + DT2overM*(CF_t.x + SF_t.x);
            // SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + (PARTPOS_T[i].y - PARTPOS_TM1[i].y)*alpha + DT2overM*(CF_t.y + SF_t.y);
            // SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + (PARTPOS_T[i].z - PARTPOS_TM1[i].z)*alpha + DT2overM*(CF_t.z + SF_t.z);

            // velocity -> momentum
            PARTMOM_T[i].x = Mi*(PARTVEL[i].x - cyclotronFreq*PARTPOS_T[i].y)*alpha;
            PARTMOM_T[i].y = Mi*(PARTVEL[i].y + cyclotronFreq*PARTPOS_T[i].x)*alpha;
            PARTMOM_T[i].z = Mi*PARTVEL[i].z*alpha;


            // half step on momentum
            // PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*(CF_t.x + SF_t.x); //forces from actual positions
            // PARTMOM_TP05[i].y = PARTMOM_T[i].y + DTover2*(CF_t.y + SF_t.y);
            // PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);
            PARTMOM_Tilday = PARTMOM_T[i].y + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_T[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
            PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*((CF_t.x + SF_t.x) + cyclotronFreq*(PARTMOM_Tilday - Mi*cyclotronFreq*PARTPOS_T[i].x)); //forces from actual positions
            PARTMOM_TP05[i].y = PARTMOM_Tilday + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP05[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
            PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);

            // // full step on positions
            // SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*overMi*PARTMOM_TP05[i].x;
            // SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + DT*overMi*PARTMOM_TP05[i].y;
            // SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;
            PARTPOS_Tilday = PARTPOS_T[i].y + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_T[i].x);
            SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*(overMi*PARTMOM_TP05[i].x + cyclotronFreq*PARTPOS_Tilday);
            SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_Tilday + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_TP1[i].x);
            SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;
          }

          //print forces
          // printf("\n Intermolecular vs Lorentz : F=%.4e \t FL=%.4e \n", sumIntermolecularForces/NPART, sumLorentzForces/NPART);
          // sumIntermolecularForces = 0;
          // sumLorentzForces = 0;


          for (i=0; i<NPART; i++) {
            Mi = M[INDX[i]];
            overMi = 1./Mi;
            cyclotronFreq = Q[INDX[i]]*B0*overMi*0.5;
            DToverMi = DT*overMi;
            DTover2 = DT*0.5;
            DTover4 = DT*0.25;
            // Recalculate forces w.r. to new positions
            if (POT == 'J') {

                CF_t = CoreForce_Jac(PARTPOS_TP1, SHELLPOS_TP1, i);
                SF_t = ShellForce_Jac(SHELLPOS_TP1, PARTPOS_TP1, i);

            } else if (POT == 'C') {

                CF_t = CoreForce_Cicc(PARTPOS_TP1, SHELLPOS_TP1, i);
                SF_t = ShellForce_Cicc(SHELLPOS_TP1, PARTPOS_TP1, i);
            }

            if (DEBUG_FLAG && _D_TOT_FORCES) {

                CFtot.x += (CF_t.x);
                CFtot.y += (CF_t.y);
                CFtot.z += (CF_t.z);

                SFtot.x += (SF_t.x);
                SFtot.y += (SF_t.y);
                SFtot.z += (SF_t.z);

                Ftot.x += (CF_t.x  + SF_t.x);
                Ftot.y += (CF_t.y  + SF_t.y);
                Ftot.z += (CF_t.z  + SF_t.z);
            }

            // final step on momentum
            // PARTMOM_TP1[i].x = PARTMOM_TP05[i].x + DTover2*(CF_t.x + SF_t.x);
            // PARTMOM_TP1[i].y = PARTMOM_TP05[i].y + DTover2*(CF_t.y + SF_t.y);
            // PARTMOM_TP1[i].z = PARTMOM_TP05[i].z + DTover2*(CF_t.z + SF_t.z);
            PARTMOM_Tilday = PARTMOM_TP05[i].y + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP05[i].x + Mi*cyclotronFreq*PARTPOS_TP1[i].y));
            PARTMOM_TP1[i].x = PARTMOM_TP05[i].x + DTover2*((CF_t.x + SF_t.x) + cyclotronFreq*(PARTMOM_Tilday - Mi*cyclotronFreq*PARTPOS_TP1[i].x));
            PARTMOM_TP1[i].y = PARTMOM_Tilday + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP1[i].x + Mi*cyclotronFreq*PARTPOS_TP1[i].y));
            PARTMOM_TP1[i].z = PARTMOM_TP05[i].z + DTover2*(CF_t.z + SF_t.z);

            // momentum -> velocity
            SHELLVEL[i].x = PARTVEL[i].x = PARTMOM_TP1[i].x*overMi + cyclotronFreq*PARTPOS_TP1[i].y;
            SHELLVEL[i].y = PARTVEL[i].y = PARTMOM_TP1[i].y*overMi - cyclotronFreq*PARTPOS_TP1[i].x;
            SHELLVEL[i].z = PARTVEL[i].z = PARTMOM_TP1[i].z*overMi;


        }

        if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

        for (i=0; i<NPART; i++) {

            //PARTVEL[i] = Velocity(PARTPOS_TM1[i], PARTPOS_TP1[i]);
            //SHELLVEL[i] = Velocity(SHELLPOS_TM1[i], SHELLPOS_TP1[i]);

            PARTPOS_TM1[i] = PARTPOS_T[i];
            PARTPOS_T[i] = PARTPOS_TP1[i];
            SHELLPOS_TM1[i] = SHELLPOS_T[i];
            SHELLPOS_T[i] = SHELLPOS_TP1[i];
        }

        if ((t+1) % IANFILE == 0) Analyse(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, therm);
        if ((t+1) % IPS == 0)  Write_PSConfig(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
        if ((t+1) % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
        if ((t+1) % IGOFR == 0) Write_GofR(t+1, PARTPOS_T);
        if ((t+1) % ICHECK == 0) Checkpoint(t+1, PARTPOS_T, SHELLPOS_T, PARTVEL, SHELLVEL);

        t_end = clock();

        Write_Elapsed_timeperstep(t+1, (double)(t_end - t_start)/CLOCKS_PER_SEC);

        if (DEBUG_FLAG && _D_PSCONFIG) {

            printf("\nPHASE-SPACE CONFIGURATION AT TIME t = %d:\n", t);
            for (i=0; i<NPART; i++) {

                printf("%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", PARTPOS_TM1[i].x, PARTPOS_TM1[i].y, PARTPOS_TM1[i].z, PARTVEL[i].x, PARTVEL[i].y, PARTVEL[i].z);
            }
        }
    }
}

void Block_MD_St_Ew(void){

    int t, i, t0 = 0;
    double DT2 = (DT*DT), DT2overM, alpha;

//    struct point Ftot = {0};
//    struct point CFtot = {0}, SFtot = {0};

    char therm;

    clock_t t_start, t_end;

    if (PRECONFIG_FLAG == 1) {

        FirstStep_St();
        t0 = 1;
    }

    for (t=t0; t<NTIMESTEPS; t++) {

        t_start = clock();

        alpha = 1.;
        therm = ' ';

//        if (DEBUG_FLAG && _D_TOT_FORCES) CFtot.x = CFtot.y = CFtot.z = SFtot.x = SFtot.y = SFtot.z = Ftot.x = Ftot.y = Ftot.z = 0;

        if (THERMOSTAT == 'T' && fabs(ITEMP - TEMP) > (_TEMP_TOL*TEMP)){

            therm = '*';
            alpha = Thermostat(PARTPOS_TM1, PARTPOS_T, TEMP);

            if (DEBUG_FLAG && _D_THERMOSTAT) {

                printf("|∆T| = %.4e\t|∆T|_M = %.4e\talpha = %.4e\n", fabs(TEMP - ITEMP), _TEMP_TOL*TEMP, alpha);
            }
        }

        Forces(PARTPOS_T, SHELLPOS_T);

        for (i=0; i<NPART; i++) {

            DT2overM = DT2/M[INDX[i]];

//            if (POT == 'J') {
//
//                CF_t = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
//                SF_t = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);
//
//            } else if (POT == 'C') {
//
//                CF_t = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
//                SF_t = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
//            }
//
//            if (DEBUG_FLAG && _D_TOT_FORCES) {
//
//                CFtot.x += (CF_t.x);
//                CFtot.y += (CF_t.y);
//                CFtot.z += (CF_t.z);
//
//                SFtot.x += (SF_t.x);
//                SFtot.y += (SF_t.y);
//                SFtot.z += (SF_t.z);
//
//                Ftot.x += (CF_t.x  + SF_t.x);
//                Ftot.y += (CF_t.y  + SF_t.y);
//                Ftot.z += (CF_t.z  + SF_t.z);
//            }

            SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + (PARTPOS_T[i].x - PARTPOS_TM1[i].x)*alpha + DT2overM*(CF[i].x + SF[i].x);
            SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + (PARTPOS_T[i].y - PARTPOS_TM1[i].y)*alpha + DT2overM*(CF[i].y + SF[i].y);
            SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + (PARTPOS_T[i].z - PARTPOS_TM1[i].z)*alpha + DT2overM*(CF[i].z + SF[i].z);
        }

//        if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

        for (i=0; i<NPART; i++) {

            PARTVEL[i] = Velocity(PARTPOS_TM1[i], PARTPOS_TP1[i]);
            SHELLVEL[i] = Velocity(SHELLPOS_TM1[i], SHELLPOS_TP1[i]);

            PARTPOS_TM1[i] = PARTPOS_T[i];
            PARTPOS_T[i] = PARTPOS_TP1[i];
            SHELLPOS_TM1[i] = SHELLPOS_T[i];
            SHELLPOS_T[i] = SHELLPOS_TP1[i];
        }

        if ((t+1) % IANFILE == 0) Analyse(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, therm);
        if ((t+1) % IPS == 0)  Write_PSConfig(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
        if ((t+1) % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
        if ((t+1) % IGOFR == 0) Write_GofR(t+1, PARTPOS_T);
        if ((t+1) % ICHECK == 0) Checkpoint(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTPOS_T, SHELLPOS_T);

        t_end = clock();

        Write_Elapsed_timeperstep(t+1, (double)(t_end - t_start)/CLOCKS_PER_SEC);

        if (DEBUG_FLAG && _D_PSCONFIG) {

            printf("\nPHASE-SPACE CONFIGURATION AT TIME t = %d:\n", t);
            for (i=0; i<NPART; i++) {

                printf("%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", PARTPOS_TM1[i].x, PARTPOS_TM1[i].y, PARTPOS_TM1[i].z, PARTVEL[i].x, PARTVEL[i].y, PARTVEL[i].z);
            }
        }
    }
}

void Block_MD_Pol(void){

    int t, i, t0 = 0;
    double DT2 = (DT*DT), DT2overM, DTover2, DTover4, overMi, Mi, alpha;
    double cyclotronFreq, PARTMOM_Tilday, PARTPOS_Tilday;
    struct point CF_t, SF_t;

    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};

    char therm;

    clock_t t_start, t_end;

    if (PRECONFIG_FLAG == 1) {

        FirstStep_Pol();
        t0 = 1;
    }
    // before the first step shell vel has to found from shellpos
    for (i=0; i<NPART; i++) {

        SHELLVEL[i] = Velocity(SHELLPOS_TM1[i], SHELLPOS_T[i]);
        SHELLVEL[i].x = SHELLVEL[i].x*2.;
        SHELLVEL[i].y = SHELLVEL[i].y*2.;
        SHELLVEL[i].z = SHELLVEL[i].z*2.;
        SHELLACC_TM1[i].x = 0;
        SHELLACC_TM1[i].y = 0;
        SHELLACC_TM1[i].z = 0;
        SHELLACC_T[i].x = 0;
        SHELLACC_T[i].y = 0;
        SHELLACC_T[i].z = 0;
        SHELLACC_TP1[i].x = 0;
        SHELLACC_TP1[i].y = 0;
        SHELLACC_TP1[i].z = 0;


    }

    for (t=t0; t<NTIMESTEPS; t++) {

        t_start = clock();

        therm = ' ';
        alpha = 1.;

        if (DEBUG_FLAG && _D_TOT_FORCES) CFtot.x = CFtot.y = CFtot.z = SFtot.x = SFtot.y = SFtot.z = Ftot.x = Ftot.y = Ftot.z = 0;

        if (THERMOSTAT == 'T' && fabs(ITEMP - TEMP) > (_TEMP_TOL*TEMP)){

            therm = '*';
            alpha = Thermostat(PARTPOS_TM1, PARTPOS_T, TEMP);

            if (DEBUG_FLAG && _D_THERMOSTAT) {

                printf("\nThermostatting: |∆T| = %.4e\t|∆T|_M = %.4e\talpha = %.4e\n\n", fabs(TEMP - ITEMP), _TEMP_TOL*TEMP, alpha);
            }
        }

        for (i=0; i<NPART; i++) {

            DT2overM = DT2/M[INDX[i]];
            Mi = M[INDX[i]];
            overMi = 1./Mi;
            cyclotronFreq = Q[INDX[i]]*B0*overMi*0.5;
            DTover2 = DT*0.5;
            DTover4 = DT*0.25;

            if (POT == 'J') {

                CF_t = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
                SF_t = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);

            } else if (POT == 'C'){

                CF_t = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
                SF_t = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
            }

            if (DEBUG_FLAG && _D_TOT_FORCES) {

                CFtot.x += (CF_t.x);
                CFtot.y += (CF_t.y);
                CFtot.z += (CF_t.z);

                SFtot.x += (SF_t.x);
                SFtot.y += (SF_t.y);
                SFtot.z += (SF_t.z);

                Ftot.x += (CF_t.x  + SF_t.x);
                Ftot.y += (CF_t.y  + SF_t.y);
                Ftot.z += (CF_t.z  + SF_t.z);
            }

            // PARTMOM_T[i].x = Mi*PARTVEL[i].x*alpha;
            // PARTMOM_T[i].y = Mi*PARTVEL[i].y*alpha;
            // PARTMOM_T[i].z = Mi*PARTVEL[i].z*alpha;
            //
            // PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*(CF_t.x + SF_t.x);
            // PARTMOM_TP05[i].y = PARTMOM_T[i].y + DTover2*(CF_t.y + SF_t.y);
            // PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);
            //
            // PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*overMi*PARTMOM_TP05[i].x;
            // PARTPOS_TP1[i].y = PARTPOS_T[i].y + DT*overMi*PARTMOM_TP05[i].y;
            // PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;
            PARTMOM_T[i].x = Mi*(PARTVEL[i].x - cyclotronFreq*PARTPOS_T[i].y)*alpha;
            PARTMOM_T[i].y = Mi*(PARTVEL[i].y + cyclotronFreq*PARTPOS_T[i].x)*alpha;
            PARTMOM_T[i].z = Mi*PARTVEL[i].z*alpha;

            PARTMOM_Tilday = PARTMOM_T[i].y + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_T[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
            PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*((CF_t.x + SF_t.x) + cyclotronFreq*(PARTMOM_Tilday - Mi*cyclotronFreq*PARTPOS_T[i].x)); //forces from actual positions
            PARTMOM_TP05[i].y = PARTMOM_Tilday + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP05[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
            PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);

            PARTPOS_Tilday = PARTPOS_T[i].y + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_T[i].x);
            PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*(overMi*PARTMOM_TP05[i].x + cyclotronFreq*PARTPOS_Tilday);
            PARTPOS_TP1[i].y = PARTPOS_Tilday + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_TP1[i].x);
            PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;

            // PARTPOS_TP1[i].x = PARTPOS_T[i].x + (PARTPOS_T[i].x - PARTPOS_TM1[i].x)*alpha + DT2overM*CF_t.x;
            // PARTPOS_TP1[i].y = PARTPOS_T[i].y + (PARTPOS_T[i].y - PARTPOS_TM1[i].y)*alpha + DT2overM*CF_t.y;
            // PARTPOS_TP1[i].z = PARTPOS_T[i].z + (PARTPOS_T[i].z - PARTPOS_TM1[i].z)*alpha + DT2overM*CF_t.z;
            //
            if (B0 == 0) {
                SHELLPOS_TP1[i].x = SHELLPOS_T[i].x + (SHELLPOS_T[i].x - SHELLPOS_TM1[i].x)*alpha;
                SHELLPOS_TP1[i].y = SHELLPOS_T[i].y + (SHELLPOS_T[i].y - SHELLPOS_TM1[i].y)*alpha;
                SHELLPOS_TP1[i].z = SHELLPOS_T[i].z + (SHELLPOS_T[i].z - SHELLPOS_TM1[i].z)*alpha;
            }
            else {
                SHELLPOS_TP1[i].x = SHELLPOS_T[i].x + SHELLVEL[i].x*DT*alpha; // alpha?
                SHELLPOS_TP1[i].y = SHELLPOS_T[i].y + SHELLVEL[i].y*DT*alpha;
                SHELLPOS_TP1[i].z = SHELLPOS_T[i].z + SHELLVEL[i].z*DT*alpha;


                if (t < 1){
                    SHELLVEL_TP1[i].x = SHELLVEL[i].x*alpha;
                    SHELLVEL_TP1[i].y = SHELLVEL[i].y*alpha;
                    SHELLVEL_TP1[i].z = SHELLVEL[i].z*alpha;
                }
                else {
                    SHELLVEL_TP1[i].x = SHELLVEL[i].x*alpha - 0.5*DT*SHELLACC_TM1[i].x;
                    SHELLVEL_TP1[i].y = SHELLVEL[i].y*alpha - 0.5*DT*SHELLACC_TM1[i].y;
                    SHELLVEL_TP1[i].z = SHELLVEL[i].z*alpha - 0.5*DT*SHELLACC_TM1[i].z;
                }
            }
        }

        if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

        if (SRMODE == 'S') {

//            ML_SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);
            if (B0 == 0){
                SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);
            }
            else{
                BSHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, SHELLVEL, SHELLVEL_TP1, 1, 0, t);
            }

        } else if (SRMODE == 'D') {

            SteepestDescent(SHELLPOS_TP1, PARTPOS_TP1);

        } else if (SRMODE == 'C') {

            ConjugateGradient(SHELLPOS_TP1, PARTPOS_TP1);
        }

        for (i=0; i<NPART; i++) {

            DT2overM = DT2/M[INDX[i]];
            Mi = M[INDX[i]];
            overMi = 1./Mi;
            cyclotronFreq = Q[INDX[i]]*B0*overMi*0.5;
            DTover2 = DT*0.5;
            DTover4 = DT*0.25;

            if (POT == 'J') {

                CF_t = CoreForce_Jac(PARTPOS_TP1, SHELLPOS_TP1, i);
                SF_t = ShellForce_Jac(SHELLPOS_TP1, PARTPOS_TP1, i);

            } else if (POT == 'C'){

                CF_t = CoreForce_Cicc(PARTPOS_TP1, SHELLPOS_TP1, i);
                SF_t = ShellForce_Cicc(SHELLPOS_TP1, PARTPOS_TP1, i);
            }

            if (DEBUG_FLAG && _D_TOT_FORCES) {

                CFtot.x += (CF_t.x);
                CFtot.y += (CF_t.y);
                CFtot.z += (CF_t.z);

                SFtot.x += (SF_t.x);
                SFtot.y += (SF_t.y);
                SFtot.z += (SF_t.z);

                Ftot.x += (CF_t.x  + SF_t.x);
                Ftot.y += (CF_t.y  + SF_t.y);
                Ftot.z += (CF_t.z  + SF_t.z);
            }

            PARTMOM_Tilday = PARTMOM_TP05[i].y + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP05[i].x + Mi*cyclotronFreq*PARTPOS_TP1[i].y));
            PARTMOM_TP1[i].x = PARTMOM_TP05[i].x + DTover2*((CF_t.x + SF_t.x) + cyclotronFreq*(PARTMOM_Tilday - Mi*cyclotronFreq*PARTPOS_TP1[i].x));
            PARTMOM_TP1[i].y = PARTMOM_Tilday + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP1[i].x + Mi*cyclotronFreq*PARTPOS_TP1[i].y));
            PARTMOM_TP1[i].z = PARTMOM_TP05[i].z + DTover2*(CF_t.z + SF_t.z);

            PARTVEL[i].x = PARTMOM_TP1[i].x*overMi + cyclotronFreq*PARTPOS_TP1[i].y;
            PARTVEL[i].y = PARTMOM_TP1[i].y*overMi - cyclotronFreq*PARTPOS_TP1[i].x;
            PARTVEL[i].z = PARTMOM_TP1[i].z*overMi;

        }

        for (i=0; i<NPART; i++) {

            //PARTVEL[i] = Velocity(PARTPOS_TM1[i], PARTPOS_TP1[i]);
            //SHELLVEL[i] = Velocity(SHELLPOS_TM1[i], SHELLPOS_TP1[i]);
            SHELLVEL[i] = SHELLVEL_TP1[i];

            SHELLACC_TM1[i] = SHELLACC_T[i];
            SHELLACC_T[i] = SHELLACC_TP1[i];
            PARTPOS_TM1[i] = PARTPOS_T[i];
            PARTPOS_T[i] = PARTPOS_TP1[i];
            SHELLPOS_TM1[i] = SHELLPOS_T[i];
            SHELLPOS_T[i] = SHELLPOS_TP1[i];
        }

        if (DEBUG_FLAG && _D_PSCONFIG) {

            printf("\nPHASE-SPACE CONFIGURATION AT TIME t = %d:\n", t+1);
            for (i=0; i<NPART; i++) {

                printf("%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", PARTPOS_T[i].x, PARTPOS_T[i].y, PARTPOS_T[i].z, PARTVEL[i].x, PARTVEL[i].y, PARTVEL[i].z);
            }
        }

        if ((t+1) % IANFILE == 0) Analyse(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, therm);
        if ((t+1) % IPS == 0)  Write_PSConfig(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
        if ((t+1) % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
        if ((t+1) % IGOFR == 0) Write_GofR(t+1, PARTPOS_T);
        if ((t+1) % ICHECK == 0) Checkpoint(t+1, PARTPOS_T, SHELLPOS_TM1, PARTVEL, SHELLPOS_T);

        t_end = clock();

        Write_Elapsed_timeperstep(t+1, (double)(t_end - t_start)/CLOCKS_PER_SEC);
    }
}

void Block_MD_Pol_Ew(void){

}
