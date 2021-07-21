//
//  FirstStep.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "FirstStep.h"

void FirstStep_St(void){

    int i;
    double DT2 = (DT*DT), DToverMi, DTover2, DTover4, overMi, Mi;
    double cyclotronFreq = 0, PARTMOM_Tilday, PARTPOS_Tilday;
    struct point vi, CF_t, SF_t;

    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};

    clock_t t_start, t_end;

    Analyse(0, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    Write_GofR(0, PARTPOS_T);

    t_start = clock();

    for (i=0; i<NPART; i++) {

        Mi = M[INDX[i]];
        overMi = 1./Mi;
        DToverMi = DT*overMi;
        DTover2 = DT*0.5;
        DTover4 = DT*0.25;

        if (POT == 'J') {

            CF_t = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
            SF_t = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);

        } else if (POT == 'C') {

            CF_t = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
            SF_t = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
        } else if (POT == 'W') {

            CF_t = CoreForce_WCA(PARTPOS_T, i);
            SF_t.x = 0;
            SF_t.y = 0;
            SF_t.z = 0;
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

        if (DEBUG_FLAG && _D_FORCES) printf("CF_t[%d] = (%.4e, %.4e, %.4e)\nSF_t[%d] = (%.4e, %.4e, %.4e)\n", i, CF_t.x, CF_t.y, CF_t.z, i, SF_t.x, SF_t.y, SF_t.z);

        // SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + PARTVEL[i].x*DT + .5*DT2overM*(CF_t.x + SF_t.x);
        // SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + PARTVEL[i].y*DT + .5*DT2overM*(CF_t.y + SF_t.y);
        // SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + PARTVEL[i].z*DT + .5*DT2overM*(CF_t.z + SF_t.z);

        // velocity -> momentum
        PARTMOM_T[i].x = Mi*PARTVEL[i].x;
        PARTMOM_T[i].y = Mi*PARTVEL[i].y;
        PARTMOM_T[i].z = Mi*PARTVEL[i].z;

        // half step on momentum
        // PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*(CF_t.x + SF_t.x);
        // PARTMOM_TP05[i].y = PARTMOM_T[i].y + DTover2*(CF_t.y + SF_t.y);
        // PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);
        PARTMOM_Tilday = PARTMOM_T[i].y + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_T[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
        PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*((CF_t.x + SF_t.x) + cyclotronFreq*(PARTMOM_Tilday - Mi*cyclotronFreq*PARTPOS_T[i].x)); //forces from actual positions
        PARTMOM_TP05[i].y = PARTMOM_Tilday + DTover4*((CF_t.y + SF_t.y) - cyclotronFreq*(PARTMOM_TP05[i].x + Mi*cyclotronFreq*PARTPOS_T[i].y));
        PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);

        // full step on positions
        // SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*overMi*PARTMOM_TP05[i].x;
        // SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + DT*overMi*PARTMOM_TP05[i].y;
        // SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;
        PARTPOS_Tilday = PARTPOS_T[i].y + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_T[i].x);
        SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*(overMi*PARTMOM_TP05[i].x + cyclotronFreq*PARTPOS_Tilday);
        SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_Tilday + DTover2*(overMi*PARTMOM_TP05[i].y - cyclotronFreq*PARTPOS_TP1[i].x);
        SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;
    }
    for (i=0; i<NPART; i++) {

        Mi = M[INDX[i]];
        overMi = 1./Mi;
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
        SHELLVEL[i].x = PARTVEL[i].x = PARTMOM_TP1[i].x*overMi;
        SHELLVEL[i].y = PARTVEL[i].y = PARTMOM_TP1[i].y*overMi;
        SHELLVEL[i].z = PARTVEL[i].z = PARTMOM_TP1[i].z*overMi;

    }

    if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

    for (i=0; i<NPART; i++) {

        // vi = Velocity(PARTPOS_T[i], PARTPOS_TP1[i]);
        //
        // PARTVEL[i].x = 2.*vi.x;
        // PARTVEL[i].y = 2.*vi.y;
        // PARTVEL[i].z = 2.*vi.z;

        PARTPOS_TM1[i] = PARTPOS_T[i];
        PARTPOS_T[i] = PARTPOS_TP1[i];
        SHELLPOS_TM1[i] = SHELLPOS_T[i];
        SHELLPOS_T[i] = SHELLPOS_TP1[i];
    }

    if (1 % IANFILE == 0) Analyse(1, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    if (1 % IPS == 0)  Write_PSConfig(1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
    if (1 % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    if (1 % IGOFR == 0) Write_GofR(1, PARTPOS_T);
    Checkpoint(0, PARTPOS_T, SHELLPOS_T, SHELLPOS_TM1, PARTVEL, SHELLVEL);

    t_end = clock();

    Write_Elapsed_timeperstep(1, (double)(t_end - t_start)/CLOCKS_PER_SEC);
}

void FirstStep_St_Ew(void){

    int i;
    double DT2 = (DT*DT), DT2overM;
    struct point vi;

    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};

    clock_t t_start, t_end;

    Analyse(0, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    Write_GofR(0, PARTPOS_T);

    t_start = clock();

    Forces(PARTPOS_T, SHELLPOS_T);

    for (i=0; i<NPART; i++) {

        DT2overM = DT2/M[INDX[i]];

//        if (POT == 'J') {
//
//            CF_t = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
//            SF_t = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);
//
//        } else if (POT == 'C') {
//
//            CF_t = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
//            SF_t = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
//        }
//
//        if (DEBUG_FLAG && _D_TOT_FORCES) {
//
//            CFtot.x += (CF_t.x);
//            CFtot.y += (CF_t.y);
//            CFtot.z += (CF_t.z);
//
//            SFtot.x += (SF_t.x);
//            SFtot.y += (SF_t.y);
//            SFtot.z += (SF_t.z);
//
//            Ftot.x += (CF_t.x  + SF_t.x);
//            Ftot.y += (CF_t.y  + SF_t.y);
//            Ftot.z += (CF_t.z  + SF_t.z);
//        }

//        if (DEBUG_FLAG && _D_FORCES) printf("CF_t[%d] = (%.4e, %.4e, %.4e)\nSF_t[%d] = (%.4e, %.4e, %.4e)\n", i, CF_t.x, CF_t.y, CF_t.z, i, SF_t.x, SF_t.y, SF_t.z);

        SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x = PARTPOS_T[i].x + PARTVEL[i].x*DT + .5*DT2overM*(CF[i].x + SF[i].x);
        SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y = PARTPOS_T[i].y + PARTVEL[i].y*DT + .5*DT2overM*(CF[i].y + SF[i].y);
        SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z = PARTPOS_T[i].z + PARTVEL[i].z*DT + .5*DT2overM*(CF[i].z + SF[i].z);
    }

    if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

    for (i=0; i<NPART; i++) {

        vi = Velocity(PARTPOS_T[i], PARTPOS_TP1[i]);

        PARTVEL[i].x = 2.*vi.x;
        PARTVEL[i].y = 2.*vi.y;
        PARTVEL[i].z = 2.*vi.z;

        PARTPOS_TM1[i] = PARTPOS_T[i];
        PARTPOS_T[i] = PARTPOS_TP1[i];
        SHELLPOS_TM1[i] = SHELLPOS_T[i];
        SHELLPOS_T[i] = SHELLPOS_TP1[i];
    }

    if (1 % IANFILE == 0) Analyse(1, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    if (1 % IPS == 0)  Write_PSConfig(1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
    if (1 % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    if (1 % IGOFR == 0) Write_GofR(1, PARTPOS_T);
    Checkpoint(0, PARTPOS_T, SHELLPOS_T, SHELLPOS_TM1, PARTVEL, SHELLVEL);

    t_end = clock();

    Write_Elapsed_timeperstep(1, (double)(t_end - t_start)/CLOCKS_PER_SEC);
}

void FirstStep_Pol(void){

    int i, indx_i;
    double DT2 = DT*DT, DT2overM, DT2overMU, DTover2, Mi, overMi;
    struct point vi, CF_t, SF_t, SF_P, SF_tp1;
    struct point velocityHalfStep;

    struct point Ftot = {0};
    struct point CFtot = {0}, SFtot = {0};

    clock_t t_start, t_end;

    Analyse(0, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    Write_GofR(0, PARTPOS_T);

    t_start = clock();
    //Conjugate Gradient to find SHELLPOS at t=0 from Lattice
    ConjugateGradient(SHELLPOS_T, PARTPOS_T);
    // Initialize particles and shells velocities to 0
    for (i=0; i<NPART; i++) {
      PARTVEL[i].x = 0;
      PARTVEL[i].y = 0;
      PARTVEL[i].z = 0;
      SHELLVEL[i].x = 0;
      SHELLVEL[i].y = 0;
      SHELLVEL[i].z = 0;
    }


    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];
        DT2overM = DT2/M[indx_i];
        DT2overMU = DT2/MU[indx_i];
        Mi = M[INDX[i]];
        overMi = 1./Mi;
        DTover2 = DT*0.5;

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

        if (DEBUG_FLAG && _D_FORCES) printf("CF_t[%d] = (%.4e, %.4e, %.4e)\nSF_t[%d] = (%.4e, %.4e, %.4e)\n", i, CF_t.x, CF_t.y, CF_t.z, i, SF_t.x, SF_t.y, SF_t.z);

        // PARTPOS_TP1[i].x = PARTPOS_T[i].x + PARTVEL[i].x*DT + .5*DT2overM*CF_t.x; //Taylor exp. to 2nd order put Velocity verlet instead
        // PARTPOS_TP1[i].y = PARTPOS_T[i].y + PARTVEL[i].y*DT + .5*DT2overM*CF_t.y;
        // PARTPOS_TP1[i].z = PARTPOS_T[i].z + PARTVEL[i].z*DT + .5*DT2overM*CF_t.z;
        // Velocity -> momentum
        PARTMOM_T[i].x = Mi*PARTVEL[i].x;
        PARTMOM_T[i].y = Mi*PARTVEL[i].y;
        PARTMOM_T[i].z = Mi*PARTVEL[i].z;
        //Half step on p
        PARTMOM_TP05[i].x = PARTMOM_T[i].x + DTover2*(CF_t.x + SF_t.x);
        PARTMOM_TP05[i].y = PARTMOM_T[i].y + DTover2*(CF_t.y + SF_t.y);
        PARTMOM_TP05[i].z = PARTMOM_T[i].z + DTover2*(CF_t.z + SF_t.z);

        PARTPOS_TP1[i].x = PARTPOS_T[i].x + DT*overMi*PARTMOM_TP05[i].x;
        PARTPOS_TP1[i].y = PARTPOS_T[i].y + DT*overMi*PARTMOM_TP05[i].y;
        PARTPOS_TP1[i].z = PARTPOS_T[i].z + DT*overMi*PARTMOM_TP05[i].z;

        if (MU[indx_i] == 0.) {

            // SHELLPOS_TP1[i].x = PARTPOS_TP1[i].x; //sp = s(0) + v(0)*dT + dT²/(2m)F(0) (last part is 0)
            // SHELLPOS_TP1[i].y = PARTPOS_TP1[i].y;
            // SHELLPOS_TP1[i].z = PARTPOS_TP1[i].z;
            //Prediction of the shells positions, Force is zero due to constraints
            SHELLPOS_TP1[i].x = SHELLPOS_T[i].x + DT*SHELLVEL[i].x;
            SHELLPOS_TP1[i].y = SHELLPOS_T[i].y + DT*SHELLVEL[i].y;
            SHELLPOS_TP1[i].z = SHELLPOS_T[i].z + DT*SHELLVEL[i].z;

        } else {

            SHELLPOS_TP1[i].x = SHELLPOS_T[i].x + PARTVEL[i].x*DT + DT2overMU*SF_t.x;
            SHELLPOS_TP1[i].y = SHELLPOS_T[i].y + PARTVEL[i].y*DT + DT2overMU*SF_t.y;
            SHELLPOS_TP1[i].z = SHELLPOS_T[i].z + PARTVEL[i].z*DT + DT2overMU*SF_t.z;
        }
    }

    if (DEBUG_FLAG && _D_TOT_FORCES) printf("\n****** CFtot = (%.4e, %.4e, %.4e) ******\n****** SFtot = (%.4e, %.4e, %.4e) ******\n****** Ftot = (%.4e, %.4e, %.4e) ******\n\n", CFtot.x, CFtot.y, CFtot.z, SFtot.x, SFtot.y, SFtot.z, Ftot.x, Ftot.y, Ftot.z);

    if (DEBUG_FLAG && _D_FORCES) {

        for (i=0; i<NPART; i++) {

            if (POT == 'J') {

                SF_P = ShellForce_Jac(SHELLPOS_TP1, PARTPOS_TP1, i);

            } else if (POT == 'C'){

                SF_P = ShellForce_Cicc(SHELLPOS_TP1, PARTPOS_TP1, i);
            }

            printf("SF_P[%d] = (%.4e, %.4e, %.4e)\n", i, SF_P.x, SF_P.y, SF_P.z);
        }
    }

    if (DEBUG_FLAG && _D_FSTEP) {

        printf("\nPOSITIONS BEFORE RELAXATION:\n");
        for (i=0; i<NPART; i++) {

            printf("r_t = (%.4e, %.4e, %.4e)\n", PARTPOS_T[i].x, PARTPOS_T[i].y, PARTPOS_T[i].z);
            printf("r_tp1 = (%.4e, %.4e, %.4e)\n", PARTPOS_TP1[i].x, PARTPOS_TP1[i].y, PARTPOS_TP1[i].z);
            printf("rho_t = (%.4e, %.4e, %.4e)\n", SHELLPOS_T[i].x, SHELLPOS_T[i].y, SHELLPOS_T[i].z);
            printf("rho_tp1 = (%.4e, %.4e, %.4e)\n\n", SHELLPOS_TP1[i].x, SHELLPOS_TP1[i].y, SHELLPOS_TP1[i].z);
        }
    }

    if (SRMODE == 'S') {

//        ML_SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);
        SHAKE(SHELLPOS_T, SHELLPOS_TP1, PARTPOS_T, PARTPOS_TP1, 1, 0);

    } else if (SRMODE == 'D') {

        SteepestDescent(SHELLPOS_TP1, PARTPOS_TP1);

    } else if (SRMODE == 'C') {

        ConjugateGradient(SHELLPOS_TP1, PARTPOS_TP1);
    }

    if (DEBUG_FLAG && _D_FSTEP) {

        printf("\nPOSITIONS AFTER FIRST STEP:\n");
        for (i=0; i<NPART; i++) {

            printf("r_t = (%.4e, %.4e, %.4e)\n", PARTPOS_T[i].x, PARTPOS_T[i].y, PARTPOS_T[i].z);
            printf("r_tp1 = (%.4e, %.4e, %.4e)\n", PARTPOS_TP1[i].x, PARTPOS_TP1[i].y, PARTPOS_TP1[i].z);
            printf("rho_t = (%.4e, %.4e, %.4e)\n", SHELLPOS_T[i].x, SHELLPOS_T[i].y, SHELLPOS_T[i].z);
            printf("rho_tp1 = (%.4e, %.4e, %.4e)\n\n", SHELLPOS_TP1[i].x, SHELLPOS_TP1[i].y, SHELLPOS_TP1[i].z);
        }
    }

    if (DEBUG_FLAG && _D_FORCES) {

        for (i=0; i<NPART; i++) {

            if (POT == 'J') {

                SF_tp1 = ShellForce_Jac(SHELLPOS_TP1, PARTPOS_TP1, i);

            } else if (POT == 'C'){

                SF_tp1 = ShellForce_Cicc(SHELLPOS_TP1, PARTPOS_TP1, i);
            }

            printf("SF_tp1[%d] = (%.4e, %.4e, %.4e)\n", i, SF_tp1.x, SF_tp1.y, SF_tp1.z);
        }
    }

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];
        DT2overM = DT2/M[indx_i];
        DT2overMU = DT2/MU[indx_i];
        Mi = M[INDX[i]];
        overMi = 1./Mi;
        DTover2 = DT*0.5;
        // Calculation of Forces with new shells positions computed with SHAKE
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

        //Full step on p
        PARTMOM_TP1[i].x = PARTMOM_TP05[i].x + DTover2*(CF_t.x + SF_t.x);
        PARTMOM_TP1[i].y = PARTMOM_TP05[i].y + DTover2*(CF_t.y + SF_t.y);
        PARTMOM_TP1[i].z = PARTMOM_TP05[i].z + DTover2*(CF_t.z + SF_t.z);
        // p -> v
        PARTVEL[i].x = PARTMOM_TP1[i].x*overMi;
        PARTVEL[i].y = PARTMOM_TP1[i].y*overMi;
        PARTVEL[i].z = PARTMOM_TP1[i].z*overMi;
    }

    for (i=0; i<NPART; i++) {

        // vi = Velocity(PARTPOS_T[i], PARTPOS_TP1[i]);
        //
        // PARTVEL[i].x = 2.*vi.x;
        // PARTVEL[i].y = 2.*vi.y;
        // PARTVEL[i].z = 2.*vi.z;

        PARTPOS_TM1[i] = PARTPOS_T[i];
        PARTPOS_T[i] = PARTPOS_TP1[i];
        SHELLPOS_TM1[i] = SHELLPOS_T[i];
        SHELLPOS_T[i] = SHELLPOS_TP1[i];
    }

    if (1 % IANFILE == 0) Analyse(1, PARTPOS_T, SHELLPOS_T, PARTVEL, ' ');
    if (1 % IPS == 0)  Write_PSConfig(1, PARTPOS_TM1, SHELLPOS_TM1, PARTVEL, SHELLVEL);
    if (1 % IVMD == 0) Write_Trajectory(PARTPOS_T, SHELLPOS_T);
    if (1 % IGOFR == 0) Write_GofR(1, PARTPOS_T);
    Checkpoint(0, PARTPOS_T, SHELLPOS_T, SHELLPOS_TM1, PARTVEL, SHELLVEL);

    t_end = clock();

    Write_Elapsed_timeperstep(1, (double)(t_end - t_start)/CLOCKS_PER_SEC);
}
