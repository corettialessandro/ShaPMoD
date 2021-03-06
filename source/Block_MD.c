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
    double DT2 = (DT*DT), DT2overM, alpha;
    struct point CF_t, SF_t;
    
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
            
            DT2overM = DT2/M[INDX[i]];
            
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

            SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + (PARTPOS_T[i].x - PARTPOS_TM1[i].x)*alpha + DT2overM*(CF_t.x + SF_t.x);
            SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + (PARTPOS_T[i].y - PARTPOS_TM1[i].y)*alpha + DT2overM*(CF_t.y + SF_t.y);
            SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + (PARTPOS_T[i].z - PARTPOS_TM1[i].z)*alpha + DT2overM*(CF_t.z + SF_t.z);
        }
        
        if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);
        
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
    double DT2 = (DT*DT), DT2overM, alpha;
    struct point CF_t, SF_t;
    
    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};
    
    char therm;
    
    clock_t t_start, t_end;
    
    if (PRECONFIG_FLAG == 1) {
        
        FirstStep_Pol();
        t0 = 1;
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

            PARTPOS_TP1[i].x = PARTPOS_T[i].x + (PARTPOS_T[i].x - PARTPOS_TM1[i].x)*alpha + DT2overM*CF_t.x;
            PARTPOS_TP1[i].y = PARTPOS_T[i].y + (PARTPOS_T[i].y - PARTPOS_TM1[i].y)*alpha + DT2overM*CF_t.y;
            PARTPOS_TP1[i].z = PARTPOS_T[i].z + (PARTPOS_T[i].z - PARTPOS_TM1[i].z)*alpha + DT2overM*CF_t.z;
            
            SHELLPOS_TP1[i].x = SHELLPOS_T[i].x + (SHELLPOS_T[i].x - SHELLPOS_TM1[i].x)*alpha;
            SHELLPOS_TP1[i].y = SHELLPOS_T[i].y + (SHELLPOS_T[i].y - SHELLPOS_TM1[i].y)*alpha;
            SHELLPOS_TP1[i].z = SHELLPOS_T[i].z + (SHELLPOS_T[i].z - SHELLPOS_TM1[i].z)*alpha;
        }
        
        if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);
        
        if (SRMODE == 'S') {
            
//            ML_SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);
            SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);

        } else if (SRMODE == 'D') {
            
            SteepestDescent(SHELLPOS_TP1, PARTPOS_TP1);
            
        } else if (SRMODE == 'C') {
            
            ConjugateGradient(SHELLPOS_TP1, PARTPOS_TP1);
        }

        for (i=0; i<NPART; i++) {
            
            PARTVEL[i] = Velocity(PARTPOS_TM1[i], PARTPOS_TP1[i]);
            SHELLVEL[i] = Velocity(SHELLPOS_TM1[i], SHELLPOS_TP1[i]);

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
        if ((t+1) % ICHECK == 0) Checkpoint(t+1, PARTPOS_TM1, SHELLPOS_TM1, PARTPOS_T, SHELLPOS_T);

        t_end = clock();
        
        Write_Elapsed_timeperstep(t+1, (double)(t_end - t_start)/CLOCKS_PER_SEC);
    }
}

void Block_MD_Pol_Ew(void){
    
}
