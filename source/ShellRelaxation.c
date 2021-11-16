//
//  ShellRelaxation.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 4/6/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "ShellRelaxation.h"

void SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount){

    int k, i, indx_i, count = 0;
    double discr = 0, kdiscr = -1.;
    double denom;
    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old;

    FILE *fp_constraints_out;
    char outputpath[_MAX_STR_LENGTH];
    sprintf(outputpath, "%sConstraints.txt", OUTPUTFOL);

    if ((fp_constraints_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    //    struct point s;

    for (k=0; k<NPART; k++) {

        if (POT == 'J') {

            Phi_old = ShellForce_Jac(rho_OLD, r_tp1, k);

        } else if (POT == 'C') {

            Phi_old = ShellForce_Cicc(rho_OLD, r_tp1, k);

        } else if (POT == 'W') {

            //Phi_old = Force_WCA(rho_OLD, k);
        }

        if (fabs(Phi_old.x) > discr) {

            discr = fabs(Phi_old.x);
            kdiscr = k + 0.1;
        }

        if (fabs(Phi_old.y) > discr) {

            discr = fabs(Phi_old.y);
            kdiscr = k + 0.2;
        }

        if (fabs(Phi_old.z) > discr) {

            discr = fabs(Phi_old.z);
            kdiscr = k + 0.3;
        }

        for (i=0; i<NPART; i++) {

            if (POT == 'J') {

                DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'C') {

                DPHIDRHO_T[k][i] = ConstTens_Cicc(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'W') {

                //DPHIDRHO_T[k][i] = ConstTens_WCA(rho_t, k, i);
            }

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) {
                printf("DPHIDRHO_T[%d][%d] =\n", k, i);
                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);
            }
        }
    }

    while (discr > LOW_TOL) { //Verifying the constraint condition

        if (VERBOSE_FLAG && _V_SHAKE){

            printf("Iteration: %d\tdiscr = (%.25e)\tkdiscr = (%.1lf)\n", count, discr, kdiscr);
        }

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_Gamma(timestep, count, GAMMA, r_tp1, rho_t, SHELLPOS_TM1);

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_S(timestep, count, rho_OLD, r_tp1);

        if (timestep % ISHAKE == 0) Write_SHAKE_output(timestep, count, discr, kdiscr);

        count++;

        if (count>_MAX_ITER) {

            //            if (GET_OUT) {
            //
            //                if (ccount == 0) printf("Stuck in a moment... Trying to get out:\n");
            //
            //                while (ccount < _MAX_ATT) {
            //
            //                    ccount++;
            //
            //                    printf("Attempt %2.d/%d: discr = %.4e\n", ccount, _MAX_ATT, discr);
            //
            //                    for (i=0; i<NPART; i++) {
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //
            //                        s.x = rho_OLD[i].x - r_tp1[i].x;
            //                        s.y = rho_OLD[i].y - r_tp1[i].y;
            //                        s.z = rho_OLD[i].z - r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        s = RNDM_Rotate(s);
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        rho_OLD[i].x = s.x + r_tp1[i].x;
            //                        rho_OLD[i].y = s.y + r_tp1[i].y;
            //                        rho_OLD[i].z = s.z + r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //                    }
            //
            //                    ML_SHAKE(rho_t, rho_OLD, r_t, r_tp1, timestep, ccount);
            //                }
            //            }

            printf("\nSHAKE.c -> SHAKE ERROR: Iteration limit exceeded. Convergence not reached!\niter = %d\tdiscr = %.10e\tkdiscr = %.1lf\n", count, discr, kdiscr);
            exit(EXIT_FAILURE);

        }else if (discr>_UP_TOL){

            printf("\nSHAKE.c -> SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
            exit(EXIT_FAILURE);
        }

        discr = 0;

        for (k=0; k<NPART; k++) { //Looping on all constraints

            //            X COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x;

            } else if (POT == 'C') {

                Phi_old.x = ShellForce_Cicc(rho_OLD, r_tp1, k).x;

            } else if (POT == 'W') {

                Phi_old.x = Force_WCA(rho_OLD, k).x;
            }

            if(fabs(Phi_old.x)>discr) {

                discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.1;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);

            for (i=0; i<NPART; i++) {

                if (POT == 'J') {

                    DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C') {

                    DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W') {

                    //DPhixDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fx;
                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);

            GAMMA[k].x = Phi_old.x/denom;
            GAMMATOT[k].x += GAMMA[k].x;

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].x = %.4e\n\n", k, GAMMA[k].x);

            for (i=0; i<NPART; i++) {

                indx_i = INDX[i];

                rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
                rho_OLD[i].y -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.y;
                rho_OLD[i].z -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.z;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NPART; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Y COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y;

            } else if (POT == 'C'){

                Phi_old.y = ShellForce_Cicc(rho_OLD, r_tp1, k).y;

            }

            if(fabs(Phi_old.y)>discr) {

                discr = fabs(Phi_old.y); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.2;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].y = %.4e\n", k, Phi_old.y);

            for (i=0; i<NPART; i++) {

                if (POT == 'J') {

                    DPhiyDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhiyDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD
                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);

                denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);

            GAMMA[k].y = Phi_old.y/denom;
            GAMMATOT[k].y += GAMMA[k].y;

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].y = %.4e\n\n", k, GAMMA[k].y);

            for (i=0; i<NPART; i++) {

                indx_i = INDX[i];

                rho_OLD[i].x -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.x;
                rho_OLD[i].y -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.y;
                rho_OLD[i].z -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.z;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NPART; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Z COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z;

            } else if (POT == 'C'){

                Phi_old.z = ShellForce_Cicc(rho_OLD, r_tp1, k).z;
            }

            if(fabs(Phi_old.z)>discr) {

                discr = fabs(Phi_old.z); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.3;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].z = %.4e\n", k, Phi_old.z);

            for (i=0; i<NPART; i++) {

                if (POT == 'J') {

                    DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhizDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);

                denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);

            GAMMA[k].z = Phi_old.z/denom;
            GAMMATOT[k].z += GAMMA[k].z;

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].z = %.4e\n\n", k, GAMMA[k].z);

            for (i=0; i<NPART; i++) {

                indx_i = INDX[i];

                rho_OLD[i].x -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.x;
                rho_OLD[i].y -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.y;
                rho_OLD[i].z -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.z;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NPART; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            if (DEBUG_FLAG && _D_CONSTR) printf("it = %d -> Phi[%d] = (%.4e, %.4e, %.4e)\n", count, k, Phi_old.x, Phi_old.y, Phi_old.z);
        } //End loop on constraints
        fprintf(fp_constraints_out, "%d \t %.10e \t %f \n", count, discr, kdiscr);



    } //End while(constraint condition)
    fprintf(fp_constraints_out, "\n");
    fflush(fp_constraints_out);
    fclose(fp_constraints_out);

    SR_ITERS = count;
    SR_DISCR = discr;

//    Uncomment for direct comparison with Conjugate Gradient method
//
//    for (i=0; i<NPART; i++) {
//
//        GAMMA[i] = ShellForce_Jac(rho_OLD, r_tp1, i);
//    }
//
//    SR_DISCR = Variance(GAMMA, NPART);

    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
}

void BSHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], struct point vrho_t[], struct point vrho_OLD[], int timestep, int ccount, int timeMD){

    int k, i, indx_i, count = 0;
    double discr = 0, kdiscr = -1.;
    double denom;
    double sumPhi = 0, meanPhi = 0;
    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old, CS_d;
    struct point DPhixDvrho_old, DPhiyDvrho_old;

    FILE *fp_constraints_out;
    char outputpath[_MAX_STR_LENGTH];
    sprintf(outputpath, "%sConstraints.txt", OUTPUTFOL);

    if ((fp_constraints_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    //    struct point s;
    for (k=0; k<NPART; k++) {

        if (POT == 'J') {

            Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x + CHI[INDX[k]]*B0*vrho_OLD[k].y;//predicted?
            Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y - CHI[INDX[k]]*B0*vrho_OLD[k].x;
            Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z; // sigma ^?
            // Phi_old.x = 0.5*Q[INDX[k]]*B0*vrho_OLD[k].y;//predicted?
            // Phi_old.y = 0.5*Q[INDX[k]]*B0*vrho_OLD[k].x;
            // Phi_old.z = 0;



        }else if (POT == 'C') {

            Phi_old = ShellForce_Cicc(rho_OLD, r_tp1, k);
        }

        if (fabs(Phi_old.x) > discr){

            discr = fabs(Phi_old.x);
            kdiscr = k + 0.1;
        }

        if (fabs(Phi_old.y) > discr) {

            discr = fabs(Phi_old.y);
            kdiscr = k + 0.2;
        }

        if (fabs(Phi_old.z) > discr) {

            discr = fabs(Phi_old.z);
            kdiscr = k + 0.3;
        }

        for (i=0; i<NPART; i++) {

            if (POT == 'J') {

                DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t) *


                DPHIDVRHO_T[k][i].fx.x = 0;
                DPHIDVRHO_T[k][i].fx.y = 0;
                DPHIDVRHO_T[k][i].fx.z = 0;

                DPHIDVRHO_T[k][i].fy.x = 0;
                DPHIDVRHO_T[k][i].fy.y = 0;
                DPHIDVRHO_T[k][i].fy.z = 0;

                DPHIDVRHO_T[k][i].fz.x = DPHIDRHO_T[k][i].fz.x;
                DPHIDVRHO_T[k][i].fz.y = DPHIDRHO_T[k][i].fz.y;
                DPHIDVRHO_T[k][i].fz.z = DPHIDRHO_T[k][i].fz.z;

                // DPHIDVRHO_T[k][i].fz.x = 0;
                // DPHIDVRHO_T[k][i].fz.y = 0;
                // DPHIDVRHO_T[k][i].fz.z = 0;

            }else if (POT == 'C') {

                DPHIDRHO_T[k][i] = ConstTens_Cicc(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)
            }



            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) {
                printf("DPHIDRHO_T[%d][%d] =\n", k, i);
                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);
            }
        }
        DPHIDVRHO_T[k][k].fx.y = CHI[INDX[k]]*B0;

        DPHIDVRHO_T[k][k].fy.x = -CHI[INDX[k]]*B0;

    }
    // printf("discr = %e\n", discr);
    for (k=0; k<NPART; k++){

        SHELLACC_TP1[k].x = 0;
        SHELLACC_TP1[k].y = 0;
        SHELLACC_TP1[k].z = 0;
    }
    //printf(" prov gamma %.4e\n", GAMMA[0].x);

    // printf("before %.20e\n", rho_OLD[0].x);
    // printf("before %.20e\n", rho_OLD[1].x);
    // printf("before %.20e\n", rho_OLD[2].x);
    // printf("before %.20e\n", rho_OLD[3].x);
    // printf("before %.20e\n", rho_OLD[4].x);
    // printf("before %.20e\n", rho_OLD[5].x);
    // printf("before %.20e\n", rho_OLD[0].y);
    // printf("before %.20e\n", rho_OLD[1].y);
    // printf("before %.20e\n", rho_OLD[2].y);
    // printf("before %.20e\n", rho_OLD[3].y);
    // printf("before %.20e\n", rho_OLD[4].y);
    // printf("before %.20e\n", rho_OLD[5].y);
    // printf("before %.20e\n", rho_OLD[0].z);
    // printf("before %.20e\n", rho_OLD[1].z);
    // printf("before %.20e\n", rho_OLD[2].z);
    // printf("before %.20e\n", rho_OLD[3].z);
    // printf("before %.20e\n", rho_OLD[4].z);
    // printf("before %.20e\n\n", rho_OLD[5].z);

    //test with old GAMMA
    for (k=0; k<NPART; k++) {
        rho_OLD[k].x -= DT*GAMMATOT[k].y*DPHIDVRHO_T[k][k].fy.x;
        rho_OLD[k].y -= DT*GAMMATOT[k].x*DPHIDVRHO_T[k][k].fx.y;
        rho_OLD[k].z -= 0;

        vrho_OLD[k].x -= (2.*GAMMATOT[k].y*DPHIDVRHO_T[k][k].fy.x);
        vrho_OLD[k].y -= (2.*GAMMATOT[k].x*DPHIDVRHO_T[k][k].fx.y);
        vrho_OLD[k].z -= 0;

        for (i=0; i<NPART; i++) {

            rho_OLD[i].x -= (DT*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.x);
            rho_OLD[i].y -= (DT*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.y);
            rho_OLD[i].z -= (DT*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.z);

            vrho_OLD[i].x -= (2.*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.x);
            vrho_OLD[i].y -= (2.*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.y);
            vrho_OLD[i].z -= (2.*GAMMATOT[k].z*DPHIDVRHO_T[k][i].fz.z);
        }
    }
    // printf("after %.20e\n", rho_OLD[0].x);
    // printf("after %.20e\n", rho_OLD[1].x);
    // printf("after %.20e\n", rho_OLD[2].x);
    // printf("after %.20e\n", rho_OLD[3].x);
    // printf("after %.20e\n", rho_OLD[4].x);
    // printf("after %.20e\n", rho_OLD[5].x);
    // printf("after %.20e\n", rho_OLD[0].y);
    // printf("after %.20e\n", rho_OLD[1].y);
    // printf("after %.20e\n", rho_OLD[2].y);
    // printf("after %.20e\n", rho_OLD[3].y);
    // printf("after %.20e\n", rho_OLD[4].y);
    // printf("after %.20e\n", rho_OLD[5].y);
    // printf("after %.20e\n", rho_OLD[0].z);
    // printf("after %.20e\n", rho_OLD[1].z);
    // printf("after %.20e\n", rho_OLD[2].z);
    // printf("after %.20e\n", rho_OLD[3].z);
    // printf("after %.20e\n", rho_OLD[4].z);
    // printf("after %.20e\n", rho_OLD[5].z);


    while (discr > LOW_TOL) { //Verifying the constraint condition

        if (VERBOSE_FLAG && _V_SHAKE){

            printf("Iteration: %d\tdiscr = (%.25e)\tkdiscr = (%.1lf)\n", count, discr, kdiscr);
        }

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_Gamma(timestep, count, GAMMA, r_tp1, rho_t, SHELLPOS_TM1);

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_S(timestep, count, rho_OLD, r_tp1);

        if (timestep % ISHAKE == 0) Write_SHAKE_output(timestep, count, discr, kdiscr);

        count++;

        if (count>_MAX_ITER) {

            //            if (GET_OUT) {
            //
            //                if (ccount == 0) printf("Stuck in a moment... Trying to get out:\n");
            //
            //                while (ccount < _MAX_ATT) {
            //
            //                    ccount++;
            //
            //                    printf("Attempt %2.d/%d: discr = %.4e\n", ccount, _MAX_ATT, discr);
            //
            //                    for (i=0; i<NPART; i++) {
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //
            //                        s.x = rho_OLD[i].x - r_tp1[i].x;
            //                        s.y = rho_OLD[i].y - r_tp1[i].y;
            //                        s.z = rho_OLD[i].z - r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        s = RNDM_Rotate(s);
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        rho_OLD[i].x = s.x + r_tp1[i].x;
            //                        rho_OLD[i].y = s.y + r_tp1[i].y;
            //                        rho_OLD[i].z = s.z + r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //                    }
            //
            //                    ML_SHAKE(rho_t, rho_OLD, r_t, r_tp1, timestep, ccount);
            //                }
            //            }

            printf("\nSHAKE.c -> SHAKE ERROR: Iteration limit exceeded. Convergence not reached!\niter = %d\tdiscr = %.10e\tkdiscr = %.1lf\n", count, discr, kdiscr);
            exit(EXIT_FAILURE);

        }else if (discr>_UP_TOL){

            printf("\nSHAKE.c -> SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
            exit(EXIT_FAILURE);
        }

        discr = 0;



        for (k=0; k<NPART; k++) { //Looping on all constraints

            //            X COMPONENT OF THE FORCE
            //denom = 0;
            if (POT == 'J') {

                Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x + CHI[INDX[k]]*B0*vrho_OLD[k].y;
                // Phi_old.x = 0.5*Q[INDX[k]]*B0*vrho_OLD[k].y;

            } else if (POT == 'C'){

                Phi_old.x = ShellForce_Cicc(rho_OLD, r_tp1, k).x;
            }

            if(fabs(Phi_old.x)>discr) {

                discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.1;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);

            //for (i=0; i<NPART; i++) {

            if (POT == 'J') {

                DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, k).fx; // TO BE COMPUTED FOR r_OLD struct?
                // DPhixDrho_old.x = 0;
                // DPhixDrho_old.y = 0;
                // DPhixDrho_old.z = 0;

                DPhixDvrho_old.y = CHI[INDX[k]]*B0;


            } else if (POT == 'C'){

                DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, k).fx; // TO BE COMPUTED FOR r_OLD
            }

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

            // if (B0 == 0){
            //     denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
            // }

                //denom += DPhixDrho_old.y*DPHIDVRHO_T[k][i].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][i].fx.y;

            denom = (DPhixDrho_old.y*DPHIDVRHO_T[k][k].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][k].fx.y*1.);

            // else {
            //     denom = (DPhixDrho_old.y*DPHIDVRHO_T[k][k].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][k].fx.y*3.);
            // }
            // denom = 0.5*DT*(DPhixDrho_old.y*DPHIDVRHO_T[k][k].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][k].fx.y*3.);
            // printf("denom.x (part %d) = %.4e \t",denom, k);

            //}

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);


            GAMMA[k].x = (SOR+SORR)*Phi_old.x/denom;
            GAMMATOT[k].x += GAMMA[k].x;
            //printf("GAMMA[%d].x = %.4e\n\n", k, GAMMA[k].x);

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].x = %.4e\n\n", k, GAMMA[k].x);

            // if (B0 == 0){
            //     for (i=0; i<NPART; i++) {
            //
            //         indx_i = INDX[i];
            //
            //         rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
            //         rho_OLD[i].y -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.y;
            //         rho_OLD[i].z -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.z;
            //     }
            //}

            // if (DEBUG_FLAG && _D_SHAKE)  {
            //
            //     for (i=0; i<NPART; i++) {
            //
            //         printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //     }
            // }

            //            Y COMPONENT OF THE FORCE
            //denom = 0;

            if (POT == 'J') {

                Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y - CHI[INDX[k]]*B0*vrho_OLD[k].x;
                // Phi_old.y = 0.5*Q[INDX[k]]*B0*vrho_OLD[k].x;

            } else if (POT == 'C'){

                Phi_old.y = ShellForce_Cicc(rho_OLD, r_tp1, k).y;
            }

            if(fabs(Phi_old.y)>discr) {

                discr = fabs(Phi_old.y); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.2;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].y = %.4e\n", k, Phi_old.y);



            if (POT == 'J') {

                DPhiyDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, k).fy; // TO BE COMPUTED FOR r_OLD
                // DPhiyDrho_old.x = 0;
                // DPhiyDrho_old.y = 0;
                // DPhiyDrho_old.z = 0;

                DPhiyDvrho_old.x = -CHI[INDX[k]]*B0;

            } else if (POT == 'C'){

                DPhiyDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, k).fy; // TO BE COMPUTED FOR r_OLD
            }

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);

            //denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);

            denom = (DPhiyDrho_old.x*DPHIDVRHO_T[k][k].fy.x*DT + DPhiyDvrho_old.x*DPHIDVRHO_T[k][k].fy.x*1.);

            // else {
            //     denom = (DPhiyDrho_old.x*DPHIDVRHO_T[k][k].fy.x*DT + DPhiyDvrho_old.x*DPHIDVRHO_T[k][k].fy.x*3.);
            // }
            //denom = 0.5*DT*(DPhiyDrho_old.x*DPHIDVRHO_T[k][k].fy.x*DT + DPhiyDvrho_old.x*DPHIDVRHO_T[k][k].fy.x*3);
            // printf("denom.y (part %d) = %.4e \t",denom, k);


            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);

            GAMMA[k].y = (SOR+SORR)*Phi_old.y/denom;
            GAMMATOT[k].y += GAMMA[k].y;
            //printf("Charge = %.4e", CHI[INDX[k]]);
            //printf("GAMMA[%d].y = %.4e\n\n", k, GAMMA[k].y);

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].y = %.4e\n\n", k, GAMMA[k].y);

            // for (i=0; i<NPART; i++) {
            //
            //     indx_i = INDX[i];
            //
            //     rho_OLD[i].x -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.x;
            //     rho_OLD[i].y -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.y;
            //     rho_OLD[i].z -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.z;
            // }

            // if (DEBUG_FLAG && _D_SHAKE)  {
            //
            //     for (i=0; i<NPART; i++) {
            //
            //         printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //     }
            // }

            //            Z COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z;
                // Phi_old.z = 0;

            } else if (POT == 'C'){

                Phi_old.z = ShellForce_Cicc(rho_OLD, r_tp1, k).z;
            }

            if(fabs(Phi_old.z)>discr) {

                discr = fabs(Phi_old.z); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.3;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].z = %.4e\n", k, Phi_old.z);

            // for (i=0; i<NPART; i++) {
            //
            //     if (POT == 'J') {
            //
            //         DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
            //
            //     } else if (POT == 'C'){
            //
            //         DPhizDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
            //     }
            //
            //     if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);
            //
            //     denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
            // }

            for (i=0; i<NPART; i++) {

                if (POT == 'J') {

                    DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
                    // DPhizDrho_old.x = 0; // TO BE COMPUTED FOR r_OLD
                    // DPhizDrho_old.y = 0;
                    // DPhizDrho_old.z = 0;

                } else if (POT == 'C'){

                    DPhizDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);

                //denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                denom += (DPhizDrho_old.x*DPHIDVRHO_T[k][i].fz.x*DT + DPhizDrho_old.y*DPHIDVRHO_T[k][i].fz.y*DT + DPhizDrho_old.z*DPHIDVRHO_T[k][i].fz.z*DT);
            }
            // printf("denom.z (part %d) = %.4e \n",denom, k);

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);

            GAMMA[k].z = Phi_old.z/denom;
            GAMMATOT[k].z += GAMMA[k].z;

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].z = %.4e\n\n", k, GAMMA[k].z);

            // for (i=0; i<NPART; i++) {
            //
            //     indx_i = INDX[i];
            //
            //     rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
            //     rho_OLD[i].y -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.y;
            //     rho_OLD[i].z -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.z;
            // }


            rho_OLD[k].x -= DT*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            rho_OLD[k].y -= DT*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            rho_OLD[k].z -= 0;

            vrho_OLD[k].x -= (2.*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x);
            vrho_OLD[k].y -= (2.*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y);
            vrho_OLD[k].z -= 0;

            // SHELLACC_TM1[k].x += GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            // SHELLACC_TM1[k].y += GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            // SHELLACC_TM1[k].z += 0;

            for (i=0; i<NPART; i++) {

                rho_OLD[i].x -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
                rho_OLD[i].y -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
                rho_OLD[i].z -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);

                vrho_OLD[i].x -= (2.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
                vrho_OLD[i].y -= (2.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
                vrho_OLD[i].z -= (2.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);

                    // SHELLACC_TM1[i].x += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x;
                    // SHELLACC_TM1[i].y += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y;
                    // SHELLACC_TM1[i].z += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z;

                    // rho_OLD[i].x -= (GAMMA[i].z*DPHIDVRHO_T[i][i].fz.x)*0.5*DT*DT);
                    // rho_OLD[i].y -= ((GAMMA[i].x*DPHIDVRHO_T[i][i].fx.y+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.y+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.y)*0.5*DT*DT);
                    // rho_OLD[i].z -= ((GAMMA[i].x*DPHIDVRHO_T[i][i].fx.z+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.z+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.z)*0.5*DT*DT);
                    //
                    // vrho_OLD[i].x -= ((GAMMA[i].x*DPHIDVRHO_T[i][i].fx.x+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.x+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.x)*1.5*DT);
                    // vrho_OLD[i].y -= ((GAMMA[i].x*DPHIDVRHO_T[i][i].fx.y+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.y+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.y)*1.5*DT);
                    // vrho_OLD[i].z -= ((GAMMA[i].x*DPHIDVRHO_T[i][i].fx.z+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.z+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.z)*1.5*DT);
                //}
            }
            // else{
            //     rho_OLD[k].x -= DT*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            //     rho_OLD[k].y -= DT*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            //     rho_OLD[k].z -= 0;
            //
            //     vrho_OLD[k].x -= (3.*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x);
            //     vrho_OLD[k].y -= (3.*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y);
            //     vrho_OLD[k].z -= 0;
            //
            //     // SHELLACC_TM1[k].x += GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            //     // SHELLACC_TM1[k].y += GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            //     // SHELLACC_TM1[k].z += 0;
            //
            //     for (i=0; i<NPART; i++) {
            //
            //         rho_OLD[i].x -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
            //         rho_OLD[i].y -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
            //         rho_OLD[i].z -= (DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);
            //
            //         vrho_OLD[i].x -= (3.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
            //         vrho_OLD[i].y -= (3.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
            //         vrho_OLD[i].z -= (3.*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);
            //     }
            // }


            // printf("Shell acc (part %d , iter %d) = (%.4e,%.4e,%.4e) \t",k,count,SHELLACC_TM1[k].x,SHELLACC_TM1[k].y,SHELLACC_TM1[k].z);
            // printf("Phi (part %d , iter %d) = (%.4e,%.4e,%.4e) \n",k,count, Phi_old.x, Phi_old.y, Phi_old.z);
            // rho_OLD[k].x -= (GAMMA[k].x*DPHIDVRHO_T[k][k].fx.x+GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x)*0.5*DT*DT;

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NPART; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            if (DEBUG_FLAG && _D_CONSTR) printf("it = %d -> Phi[%d] = (%.4e, %.4e, %.4e)\n", count, k, Phi_old.x, Phi_old.y, Phi_old.z);
        //printf("gamma = %.4e %.4e %.4e\n", GAMMA[0].x, GAMMA[0].y, GAMMA[0].z);
        // exit(0);
        //sumPhi += sqrt(Phi_old.x*Phi_old.x + Phi_old.y*Phi_old.y + Phi_old.z*Phi_old.z);




        } //End loop on constraints
        //printf("nb of iter = %d,\t discr = %e, \t discrk = %.1lf \n", count, discr,kdiscr);

        //fprintf(fp_constraints_out, "%d \t %.10e \t %f \n", count, discr, kdiscr);



    } //End while(constraint condition)
    // printf("after iteration %.20e\n", rho_OLD[0].x);
    // printf("after iteration %.20e\n", rho_OLD[1].x);
    // printf("after iteration %.20e\n", rho_OLD[2].x);
    // printf("after iteration %.20e\n", rho_OLD[3].x);
    // printf("after iteration %.20e\n", rho_OLD[4].x);
    // printf("after iteration %.20e\n", rho_OLD[5].x);
    // printf("after iteration %.20e\n", rho_OLD[0].y);
    // printf("after iteration %.20e\n", rho_OLD[1].y);
    // printf("after iteration %.20e\n", rho_OLD[2].y);
    // printf("after iteration %.20e\n", rho_OLD[3].y);
    // printf("after iteration %.20e\n", rho_OLD[4].y);
    // printf("after iteration %.20e\n", rho_OLD[5].y);
    // printf("after iteration %.20e\n", rho_OLD[0].z);
    // printf("after iteration %.20e\n", rho_OLD[1].z);
    // printf("after iteration %.20e\n", rho_OLD[2].z);
    // printf("after iteration %.20e\n", rho_OLD[3].z);
    // printf("after iteration %.20e\n", rho_OLD[4].z);
    // printf("after iteration %.20e\n", rho_OLD[5].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[0].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[1].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[2].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[3].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[4].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[5].x);
    // printf("end loop gamma %.4e\n", GAMMATOT[0].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[1].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[2].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[3].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[4].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[5].y);
    // printf("end loop gamma %.4e\n", GAMMATOT[0].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[1].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[2].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[3].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[4].z);
    // printf("end loop gamma %.4e\n", GAMMATOT[5].z);
    fprintf(fp_constraints_out, "\n");
    fflush(fp_constraints_out);
    fclose(fp_constraints_out);
    //Calculate Shell acceleration needed for the next provisional
    for (k=0; k<NPART; k++) {
        SHELLACC_TP1[k].x += GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
        SHELLACC_TP1[k].y += GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
        SHELLACC_TP1[k].z += 0;


        for (i=0; i<NPART; i++) {

            SHELLACC_TP1[i].x += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x;
            SHELLACC_TP1[i].y += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y;
            SHELLACC_TP1[i].z += GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z;
        }
        //printf("gamma = %.4e %.4e %.4e\n", GAMMA[0].x, GAMMA[0].y, GAMMA[0].z);
        //printf("Shell acc t+1 = %.4e %.4e %.4e\n", SHELLACC_TP1[0].x, SHELLACC_TP1[0].y, SHELLACC_TP1[0].z);
    }
    //printf("Shell acc t+1 = %.4e %.4e %.4e\n", SHELLACC_TP1[0].x, SHELLACC_TP1[0].y, SHELLACC_TP1[0].z);
    //printf("Acc part 1 = (%.4e, %.4e, %.4e) \n", SHELLACC_TM1[1].x, SHELLACC_TM1[1].y, SHELLACC_TM1[1].z);



    // CS_d.x = (r_tp1[1].x - rho_OLD[1].x);
    // CS_d.y = (r_tp1[1].y - rho_OLD[1].y);
    // CS_d.z = (r_tp1[1].z - rho_OLD[1].z);
    // printf("Distance part 1 = (%.4e, %.4e, %.4e) \n", CS_d.x, CS_d.y, CS_d.z);


    SR_ITERS = count;
    SR_DISCR = discr;

//    Uncomment for direct comparison with Conjugate Gradient method
//
//    for (i=0; i<NPART; i++) {
//
//        GAMMA[i] = ShellForce_Jac(rho_OLD, r_tp1, i);
//    }
//
//    SR_DISCR = Variance(GAMMA, NPART);

    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
}

void MultiSHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount){

    int k, i, j, indx_i, count = 0;
    double discr = 0, kdiscr = -1.;
    double denom, denomx, denomy, denomz, CC_r;
    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old, testForce, testdForcex, testdForcey, testdForcez, CC_d;

    FILE *fp_constraints_out;
    char outputpath[_MAX_STR_LENGTH];
    sprintf(outputpath, "%sConstraints.txt", OUTPUTFOL);

    if ((fp_constraints_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    //    struct point s;

    for (i=NATOMSPERSPEC[1]; i<NPART; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(r_t[i],i);
    }

    for (i=0; i<NATOMSPERSPEC[1]; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(rho_t[i],i);
    }

    for (k=0; k<NATOMSPERSPEC[1]; k++) {

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            if (POT == 'J') {

                DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'C') {

                DPHIDRHO_T[k][i] = ConstTens_Cicc(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'W') {

                DPHIDRHO_T[k][i] = ConstTens_WCA(rho_t, r_t, k, i);
                // DPHIDRHO_T[k][i].fx.x += 1e-8;
                // DPHIDRHO_T[k][i].fx.y += 1e-8;
                // DPHIDRHO_T[k][i].fx.z += 1e-8;
                // DPHIDRHO_T[k][i].fy.x += 1e-8;
                // DPHIDRHO_T[k][i].fy.y += 1e-8;
                // DPHIDRHO_T[k][i].fy.z += 1e-8;
                // DPHIDRHO_T[k][i].fz.x += 1e-8;
                // DPHIDRHO_T[k][i].fz.y += 1e-8;
                // DPHIDRHO_T[k][i].fz.z += 1e-8;

            } else if (POT == 'L') {

                DPHIDRHO_T[k][i] = ConstTens_LJ(rho_t, r_t, k, i);

            }
            //if (k==i)printf(" %d / %d dFxold = %.4e %.4e %.4e \n dFyold = %.4e %.4e %.4e \n dFzold = %.4e %.4e %.4e\n\n", k, i,  DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) {
                printf("DPHIDRHO_T[%d][%d] =\n", k, i);
                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);
            }
        }
    }

    for (i=NATOMSPERSPEC[1]; i<NPART; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(r_tp1[i],i);
    }

    for (i=0; i<NATOMSPERSPEC[1]; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(rho_OLD[i],i);
    }
    for (k=0; k<NATOMSPERSPEC[1]; k++) {

        if (POT == 'J') {

            Phi_old = ShellForce_Jac(rho_OLD, r_tp1, k);

        } else if (POT == 'C') {

            Phi_old = ShellForce_Cicc(rho_OLD, r_tp1, k);

        } else if (POT == 'W') {

            Phi_old = Force_WCA(rho_OLD, k);
            //printf("%.4e %.4e %.4e \n", Phi_old.x, Phi_old.y, Phi_old.z);
        } else if (POT == 'L') {

            Phi_old = Force_LJ(rho_OLD, k);
            //printf("%.4e %.4e %.4e \n", Phi_old.x, Phi_old.y, Phi_old.z);
        }

        if (fabs(Phi_old.x) > discr) {

            discr = fabs(Phi_old.x);
            kdiscr = k + 0.1;
        }

        if (fabs(Phi_old.y) > discr) {

            discr = fabs(Phi_old.y);
            kdiscr = k + 0.2;
        }

        if (fabs(Phi_old.z) > discr) {

            discr = fabs(Phi_old.z);
            kdiscr = k + 0.3;
        }
        printf("%.4e \n", Phi_old.x);
    }

    //print full SHAKE matrix w.r. to x coordinate

    for (k=0; k<NATOMSPERSPEC[1]; k++) {


        for (j=0; j<NATOMSPERSPEC[1]; j++) {
            denomx = 0.;
            denomy = 0.;
            denomz = 0.;

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'W') {

                    DPhixDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fx;
                    DPhiyDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fy;
                    DPhizDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fz;
                }else if (POT == 'L') {

                    DPhixDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fx;
                    DPhiyDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fy;
                    DPhizDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fz;
                }

                denomx += (DPHIDRHO_T[k][i].fx.x*DPHIDRHO_T[j][i].fx.x + DPHIDRHO_T[k][i].fx.y*DPHIDRHO_T[j][i].fx.y + DPHIDRHO_T[k][i].fx.z*DPHIDRHO_T[j][i].fx.z);
                // denomx += (DPhixDrho_old.x*DPHIDRHO_T[j][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[j][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[j][i].fx.z);
                // denomy += (DPhiyDrho_old.x*DPHIDRHO_T[j][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[j][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[j][i].fy.z);
                // denomz += (DPhizDrho_old.x*DPHIDRHO_T[j][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[j][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[j][i].fz.z);

            }

            printf(" %.4e \t", k, j, denomx);
            //printf("Sx(%d, %d) = %.4e \t", k, j, denomx);
            // printf("Sy(%d, %d) = %.4e \t", k, j, denomy);
            // printf("Sz(%d, %d) = %.4e \t", k, j, denomz);
            // CC_d = Distance(rho_OLD[k], rho_OLD[j]);
            // CC_r = mod(CC_d);
            // printf("dist = %.4e\n",CC_r);


        }
        printf("\n");
    }
    exit(0);

    while (discr > LOW_TOL) { //Verifying the constraint condition

        if (VERBOSE_FLAG && _V_SHAKE){

            printf("Iteration: %d\tdiscr = (%.25e)\tkdiscr = (%.1lf)\n", count, discr, kdiscr);
        }

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_Gamma(timestep, count, GAMMA, r_tp1, rho_t, SHELLPOS_TM1);

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_S(timestep, count, rho_OLD, r_tp1);

        if (timestep % ISHAKE == 0) Write_SHAKE_output(timestep, count, discr, kdiscr);

        count++;

        if (count>_MAX_ITER) {

            //            if (GET_OUT) {
            //
            //                if (ccount == 0) printf("Stuck in a moment... Trying to get out:\n");
            //
            //                while (ccount < _MAX_ATT) {
            //
            //                    ccount++;
            //
            //                    printf("Attempt %2.d/%d: discr = %.4e\n", ccount, _MAX_ATT, discr);
            //
            //                    for (i=0; i<NPART; i++) {
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //
            //                        s.x = rho_OLD[i].x - r_tp1[i].x;
            //                        s.y = rho_OLD[i].y - r_tp1[i].y;
            //                        s.z = rho_OLD[i].z - r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        s = RNDM_Rotate(s);
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        rho_OLD[i].x = s.x + r_tp1[i].x;
            //                        rho_OLD[i].y = s.y + r_tp1[i].y;
            //                        rho_OLD[i].z = s.z + r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //                    }
            //
            //                    ML_SHAKE(rho_t, rho_OLD, r_t, r_tp1, timestep, ccount);
            //                }
            //            }

            printf("\nSHAKE.c -> SHAKE ERROR: Iteration limit exceeded. Convergence not reached!\niter = %d\tdiscr = %.10e\tkdiscr = %.1lf\n", count, discr, kdiscr);
            exit(EXIT_FAILURE);

        }else if (discr>_UP_TOL){

            printf("\nSHAKE.c -> SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
            exit(EXIT_FAILURE);
        }

        discr = 0;

        for (k=0; k<NATOMSPERSPEC[1]; k++) { //Looping on all constraints

            //            X COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x;

            } else if (POT == 'C') {

                Phi_old.x = ShellForce_Cicc(rho_OLD, r_tp1, k).x;

            } else if (POT == 'W') {

                Phi_old.x = Force_WCA(rho_OLD, k).x;

            } else if (POT == 'L') {

                Phi_old.x = Force_LJ(rho_OLD, k).x;
            }

            if(fabs(Phi_old.x)>discr) {

                discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.1;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C') {

                    DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W') {
                    //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhixDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fx;
                    //if (k ==662) printf("(%f %f %f)\n",DPhixDrho_old.x,DPhixDrho_old.y,DPhixDrho_old.z);

                } else if (POT == 'L') {

                    DPhixDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fx;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

                // if (DPhixDrho_old.x != 0. && DPHIDRHO_T[k][i].fx.x == 0.) DPHIDRHO_T[k][i].fx.x =1.;
                // if (DPhixDrho_old.y != 0. && DPHIDRHO_T[k][i].fx.y == 0.) DPHIDRHO_T[k][i].fx.y =1.;
                // if (DPhixDrho_old.z != 0. && DPHIDRHO_T[k][i].fx.z == 0.) DPHIDRHO_T[k][i].fx.z =1.;
                // if (k == 42) printf("denom = %.8e\n", denom);
                // if (k == 42) printf("(%d)+= %.8e\n", i,(DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z));
                // if (DPhixDrho_old.x != 0. && DPHIDRHO_T[k][i].fx.x == 0. && DPhixDrho_old.y != 0. && DPHIDRHO_T[k][i].fx.y == 0. && DPhixDrho_old.z != 0. && DPHIDRHO_T[k][i].fx.z == 0.) {
                //     denom += fabs(DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                // } else {
                //     denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                // }


                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                //denom += (DPhixDrho_old.x*DPhixDrho_old.x + DPhixDrho_old.y*DPhixDrho_old.y + DPhixDrho_old.z*DPhixDrho_old.z);
            }
            //if (k==662) printf("%.4e\n", denom);

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);

            if (denom != 0.){
                GAMMA[k].x = (SOR)*Phi_old.x/denom;
                GAMMATOT[k].x += GAMMA[k].x;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].x = %.4e\n\n", k, GAMMA[k].x);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];

                    rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
                    rho_OLD[i].y -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.y;
                    rho_OLD[i].z -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.z;
                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].x = 0.;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Y COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y;

            } else if (POT == 'C'){

                Phi_old.y = ShellForce_Cicc(rho_OLD, r_tp1, k).y;

            } else if (POT == 'W') {

                Phi_old.y = Force_WCA(rho_OLD, k).y;

            } else if (POT == 'L') {

                Phi_old.y = Force_LJ(rho_OLD, k).y;

            }

            if(fabs(Phi_old.y)>discr) {

                discr = fabs(Phi_old.y); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.2;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].y = %.4e\n", k, Phi_old.y);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhiyDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhiyDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W') {

                  //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhiyDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fy;
                  //printf("(%f %f %f)\n",DPhiyDrho_old.x,DPhiyDrho_old.y,DPhiyDrho_old.z);

                } else if (POT == 'L') {

                    DPhiyDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fy;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);

                // if (DPhiyDrho_old.x != 0 && DPHIDRHO_T[k][i].fy.x == 0) DPHIDRHO_T[k][i].fy.x =1.;
                // if (DPhiyDrho_old.y != 0 && DPHIDRHO_T[k][i].fy.y == 0) DPHIDRHO_T[k][i].fy.y =1.;
                // if (DPhiyDrho_old.z != 0 && DPHIDRHO_T[k][i].fy.z == 0) DPHIDRHO_T[k][i].fy.z =1.;

                // if (DPhiyDrho_old.x != 0. && DPHIDRHO_T[k][i].fy.x == 0. && DPhiyDrho_old.y != 0. && DPHIDRHO_T[k][i].fy.y == 0. && DPhiyDrho_old.z != 0. && DPHIDRHO_T[k][i].fy.z == 0.) {
                //     denom += fabs(DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                // } else {
                //
                //     denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                // }

                denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                //denom += (DPhiyDrho_old.x*DPhiyDrho_old.x + DPhiyDrho_old.y*DPhiyDrho_old.y + DPhiyDrho_old.z*DPhiyDrho_old.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);

            if (denom != 0.){
                GAMMA[k].y = (SOR)*Phi_old.y/denom;
                GAMMATOT[k].y += GAMMA[k].y;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].y = %.4e\n\n", k, GAMMA[k].y);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];
                    //printf("%.4e \n", GAMMA[k].y);
                    rho_OLD[i].x -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.x; //segmentation fault
                    rho_OLD[i].y -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.y;
                    rho_OLD[i].z -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.z;

                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].y = 0.;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Z COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z;

            } else if (POT == 'C'){

                Phi_old.z = ShellForce_Cicc(rho_OLD, r_tp1, k).z;

            } else if (POT == 'W'){

                Phi_old.z = Force_WCA(rho_OLD, k).z;

            } else if (POT == 'L'){

                Phi_old.z = Force_LJ(rho_OLD, k).z;

            }

            if(fabs(Phi_old.z)>discr) {

                discr = fabs(Phi_old.z); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.3;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].z = %.4e\n", k, Phi_old.z);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhizDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W'){

                    //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhizDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fz;
                    //printf("(%f %f %f)\n",DPhizDrho_old.x,DPhizDrho_old.y,DPhizDrho_old.z);
                }  else if (POT == 'L'){

                    DPhizDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fz;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);

                // if (DPhizDrho_old.x != 0 && DPHIDRHO_T[k][i].fz.x == 0) DPHIDRHO_T[k][i].fz.x =1.;
                // if (DPhizDrho_old.y != 0 && DPHIDRHO_T[k][i].fz.y == 0) DPHIDRHO_T[k][i].fz.y =1.;
                // if (DPhizDrho_old.z != 0 && DPHIDRHO_T[k][i].fz.z == 0) DPHIDRHO_T[k][i].fz.z =1.;
                //if (k==672)printf("%.8e\n",denom);
                // if (DPhizDrho_old.x != 0. && DPHIDRHO_T[k][i].fz.x == 0. && DPhizDrho_old.y != 0. && DPHIDRHO_T[k][i].fz.y == 0. && DPhizDrho_old.z != 0. && DPHIDRHO_T[k][i].fz.z == 0.) {
                //     denom += fabs(DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                // } else {
                //     denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                // }


                denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                //denom += (DPhizDrho_old.x*DPhizDrho_old.x + DPhizDrho_old.y*DPhizDrho_old.y + DPhizDrho_old.z*DPhizDrho_old.z);
                // if (k == 672) printf("denom = %.8e\n", denom);
                // if (k == 672) printf("(%d)+= %.8e\n", i,(DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z));
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);
            //if (count == 2 && k == 448) printf("%.4e\n",denom);
            if (denom != 0.){

                GAMMA[k].z = (SOR)*Phi_old.z/denom;
                GAMMATOT[k].z += GAMMA[k].z;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].z = %.4e\n\n", k, GAMMA[k].z);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];

                    rho_OLD[i].x -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.x;
                    rho_OLD[i].y -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.y;
                    rho_OLD[i].z -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.z;

                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].z = 0.;
            }
            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            if (DEBUG_FLAG && _D_CONSTR) printf("it = %d -> Phi[%d] = (%.4e, %.4e, %.4e)\n", count, k, Phi_old.x, Phi_old.y, Phi_old.z);
            //printf("part : %d \n", k);
        } //End loop on constraints
        //if (count == 2) exit(0);
        fprintf(fp_constraints_out, "%d \t %.10e \t %f \n", count, discr, kdiscr);
        testForce = Force_WCA(rho_OLD, 171);
        testdForcex = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fx;
        testdForcey = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fy;
        testdForcez = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fz;

        printf("nb of iter = %d,\t discr = %e, \t discrk = %.1lf \n", count, discr,kdiscr);
        // printf("indx = %d, shellpos = %.4e %.4e %.4e \n, force = %.4e %.4e %.4e \n",171, rho_OLD[171].x, rho_OLD[171].y, rho_OLD[171].z, testForce.x, testForce.y, testForce.z);
        // printf("dFx = %.4e %.4e %.4e \n dFy = %.4e %.4e %.4e \n dFz = %.4e %.4e %.4e \n", testdForcex.x, testdForcex.y, testdForcex.z, testdForcey.x, testdForcey.y, testdForcey.z, testdForcez.x, testdForcez.y, testdForcez.z);
        // printf("dFxold = %.4e %.4e %.4e \n dFyold = %.4e %.4e %.4e \n dFzold = %.4e %.4e %.4e\n\n", DPHIDRHO_T[171][171].fx.x, DPHIDRHO_T[171][171].fx.y, DPHIDRHO_T[171][171].fx.z, DPHIDRHO_T[171][171].fy.x, DPHIDRHO_T[171][171].fy.y, DPHIDRHO_T[171][171].fy.z, DPHIDRHO_T[171][171].fz.x, DPHIDRHO_T[171][171].fz.y, DPHIDRHO_T[171][171].fz.z);
    } //End while(constraint condition)
    fprintf(fp_constraints_out, "\n");
    fflush(fp_constraints_out);
    fclose(fp_constraints_out);

    SR_ITERS = count;
    SR_DISCR = discr;

//    Uncomment for direct comparison with Conjugate Gradient method
//
//    for (i=0; i<NPART; i++) {
//
//        GAMMA[i] = ShellForce_Jac(rho_OLD, r_tp1, i);
//    }
//
//    SR_DISCR = Variance(GAMMA, NPART);

    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
}

void MultiWeinbachElber(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount){

    int k, i, j, indx_i, count = 0;
    double discr = 0, kdiscr = -1.;
    double denom, denomx, denomy, denomz, CC_r;
    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old, testForce, testdForcex, testdForcey, testdForcez, CC_d;

    FILE *fp_constraints_out;
    char outputpath[_MAX_STR_LENGTH];
    sprintf(outputpath, "%sConstraints.txt", OUTPUTFOL);

    if ((fp_constraints_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    //    struct point s;

    for (i=NATOMSPERSPEC[1]; i<NPART; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(r_t[i],i);
    }

    for (i=0; i<NATOMSPERSPEC[1]; i++) {
        Rem_Point_From_Cell(i);
        Add_Point_To_Cell(rho_t[i],i);
    }

    for (k=0; k<NATOMSPERSPEC[1]; k++) {

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            if (POT == 'J') {

                DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'C') {

                DPHIDRHO_T[k][i] = ConstTens_Cicc(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            } else if (POT == 'W') {

                DPHIDRHO_T[k][i] = ConstTens_WCA(rho_t, r_t, k, i);

            } else if (POT == 'L') {

                DPHIDRHO_T[k][i] = ConstTens_LJ(rho_t, r_t, k, i);

            }
            //if (k==i)printf(" %d / %d dFxold = %.4e %.4e %.4e \n dFyold = %.4e %.4e %.4e \n dFzold = %.4e %.4e %.4e\n\n", k, i,  DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) {
                printf("DPHIDRHO_T[%d][%d] =\n", k, i);
                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);
            }
        }
    }

    // Compute symmetry Shake matrix for each component
    for (k=0; k<NATOMSPERSPEC[1]; k++) {
        for (i=0; i<NATOMSPERSPEC[1]; i++) {
            for (j=0; j<NATOMSPERSPEC[1]; j++) {

                SHAKEMATRIX_X[k][i] += DPHIDRHO_T[k][j].fx.x*DPHIDRHO_T[i][j].fx.x + DPHIDRHO_T[k][j].fx.y*DPHIDRHO_T[i][j].fx.y + DPHIDRHO_T[k][j].fx.z*DPHIDRHO_T[i][j].fx.z;
                SHAKEMATRIX_Y[k][i] += DPHIDRHO_T[k][j].fy.x*DPHIDRHO_T[i][j].fy.x + DPHIDRHO_T[k][j].fy.y*DPHIDRHO_T[i][j].fy.y + DPHIDRHO_T[k][j].fy.z*DPHIDRHO_T[i][j].fy.z;
                SHAKEMATRIX_Z[k][i] += DPHIDRHO_T[k][j].fz.x*DPHIDRHO_T[i][j].fz.x + DPHIDRHO_T[k][j].fz.y*DPHIDRHO_T[i][j].fz.y + DPHIDRHO_T[k][j].fz.z*DPHIDRHO_T[i][j].fz.z;


            }
        printf("%.4e ", SHAKEMATRIX_X[k][i]);
        }
    printf("\n");
    }

    for (k=0; k<NATOMSPERSPEC[1]; k++) {

        if (POT == 'J') {

            Phi_old = ShellForce_Jac(rho_OLD, r_t, k);

        } else if (POT == 'C') {

            Phi_old = ShellForce_Cicc(rho_OLD, r_t, k);

        } else if (POT == 'W') {

            Phi_old = Force_WCA(rho_OLD, k);
            //printf("%.4e %.4e %.4e \n", Phi_old.x, Phi_old.y, Phi_old.z);
        } else if (POT == 'L') {

            Phi_old = Force_LJ(rho_OLD, k);
            //printf("%.4e %.4e %.4e \n", Phi_old.x, Phi_old.y, Phi_old.z);
        }

        PHI[k].x = Phi_old.x;
        PHI[k].y = Phi_old.y;
        PHI[k].z = Phi_old.z;

        //printf("%.4e\n", PHI[k].x);


        if (fabs(Phi_old.x) > discr) {

            discr = fabs(Phi_old.x);
            kdiscr = k + 0.1;
        }

        if (fabs(Phi_old.y) > discr) {

            discr = fabs(Phi_old.y);
            kdiscr = k + 0.2;
        }

        if (fabs(Phi_old.z) > discr) {

            discr = fabs(Phi_old.z);
            kdiscr = k + 0.3;
        }

    }

    while (discr > LOW_TOL) { //Verifying the constraint condition

        if (VERBOSE_FLAG && _V_SHAKE){

            printf("Iteration: %d\tdiscr = (%.25e)\tkdiscr = (%.1lf)\n", count, discr, kdiscr);
        }

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_Gamma(timestep, count, GAMMA, r_tp1, rho_t, SHELLPOS_TM1);

        if (_O_SHAKE && timestep % ISHAKE == 0 && count != 0) Write_S(timestep, count, rho_OLD, r_tp1);

        if (timestep % ISHAKE == 0) Write_SHAKE_output(timestep, count, discr, kdiscr);

        count++;

        if (count>_MAX_ITER) {

            //            if (GET_OUT) {
            //
            //                if (ccount == 0) printf("Stuck in a moment... Trying to get out:\n");
            //
            //                while (ccount < _MAX_ATT) {
            //
            //                    ccount++;
            //
            //                    printf("Attempt %2.d/%d: discr = %.4e\n", ccount, _MAX_ATT, discr);
            //
            //                    for (i=0; i<NPART; i++) {
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //
            //                        s.x = rho_OLD[i].x - r_tp1[i].x;
            //                        s.y = rho_OLD[i].y - r_tp1[i].y;
            //                        s.z = rho_OLD[i].z - r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("OLD: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        s = RNDM_Rotate(s);
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: s[%d] = (%.4e, %.4e, %.4e)\n", i, s.x, s.y, s.z);
            //
            //                        rho_OLD[i].x = s.x + r_tp1[i].x;
            //                        rho_OLD[i].y = s.y + r_tp1[i].y;
            //                        rho_OLD[i].z = s.z + r_tp1[i].z;
            //
            //                        if (DEBUG_FLAG && _D_STUCK) printf("NEW: rho_old[%d] = (%lf, %lf, %lf)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
            //                    }
            //
            //                    ML_SHAKE(rho_t, rho_OLD, r_t, r_tp1, timestep, ccount);
            //                }
            //            }

            printf("\nSHAKE.c -> SHAKE ERROR: Iteration limit exceeded. Convergence not reached!\niter = %d\tdiscr = %.10e\tkdiscr = %.1lf\n", count, discr, kdiscr);
            exit(EXIT_FAILURE);

        }else if (discr>_UP_TOL){

            printf("\nSHAKE.c -> SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
            exit(EXIT_FAILURE);
        }

        discr = 0;

        // SHAKEMATRIX_X[0][0] = 4.;
        // SHAKEMATRIX_X[0][1] = 1.;
        // SHAKEMATRIX_X[1][0] = 1.;
        // SHAKEMATRIX_X[1][1] = 3.;
        //
        // PHI[0].x = 1.;
        // PHI[1].x = 2.;
        //
        // GAMMA[0].x = 2.;
        // GAMMA[1].x = 1.;

        LinearConjugateGradient(SHAKEMATRIX_X, PHI, GAMMA, NATOMSPERSPEC[1], 'x');

        exit(0);

        for (k=0; k<NATOMSPERSPEC[1]; k++) { //Looping on all constraints

            //            X COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x;

            } else if (POT == 'C') {

                Phi_old.x = ShellForce_Cicc(rho_OLD, r_tp1, k).x;

            } else if (POT == 'W') {

                Phi_old.x = Force_WCA(rho_OLD, k).x;

            } else if (POT == 'L') {

                Phi_old.x = Force_LJ(rho_OLD, k).x;
            }

            if(fabs(Phi_old.x)>discr) {

                discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.1;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C') {

                    DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W') {
                    //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhixDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fx;
                    //if (k ==662) printf("(%f %f %f)\n",DPhixDrho_old.x,DPhixDrho_old.y,DPhixDrho_old.z);

                } else if (POT == 'L') {

                    DPhixDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fx;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

                // if (DPhixDrho_old.x != 0. && DPHIDRHO_T[k][i].fx.x == 0.) DPHIDRHO_T[k][i].fx.x =1.;
                // if (DPhixDrho_old.y != 0. && DPHIDRHO_T[k][i].fx.y == 0.) DPHIDRHO_T[k][i].fx.y =1.;
                // if (DPhixDrho_old.z != 0. && DPHIDRHO_T[k][i].fx.z == 0.) DPHIDRHO_T[k][i].fx.z =1.;
                // if (k == 42) printf("denom = %.8e\n", denom);
                // if (k == 42) printf("(%d)+= %.8e\n", i,(DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z));
                // if (DPhixDrho_old.x != 0. && DPHIDRHO_T[k][i].fx.x == 0. && DPhixDrho_old.y != 0. && DPHIDRHO_T[k][i].fx.y == 0. && DPhixDrho_old.z != 0. && DPHIDRHO_T[k][i].fx.z == 0.) {
                //     denom += fabs(DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                // } else {
                //     denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                // }


                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
                //denom += (DPhixDrho_old.x*DPhixDrho_old.x + DPhixDrho_old.y*DPhixDrho_old.y + DPhixDrho_old.z*DPhixDrho_old.z);
            }
            //if (k==662) printf("%.4e\n", denom);

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);

            if (denom != 0.){
                GAMMA[k].x = (SOR)*Phi_old.x/denom;
                GAMMATOT[k].x += GAMMA[k].x;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].x = %.4e\n\n", k, GAMMA[k].x);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];

                    rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
                    rho_OLD[i].y -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.y;
                    rho_OLD[i].z -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.z;
                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].x = 0.;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Y COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y;

            } else if (POT == 'C'){

                Phi_old.y = ShellForce_Cicc(rho_OLD, r_tp1, k).y;

            } else if (POT == 'W') {

                Phi_old.y = Force_WCA(rho_OLD, k).y;

            } else if (POT == 'L') {

                Phi_old.y = Force_LJ(rho_OLD, k).y;

            }

            if(fabs(Phi_old.y)>discr) {

                discr = fabs(Phi_old.y); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.2;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].y = %.4e\n", k, Phi_old.y);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhiyDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhiyDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W') {

                  //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhiyDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fy;
                  //printf("(%f %f %f)\n",DPhiyDrho_old.x,DPhiyDrho_old.y,DPhiyDrho_old.z);

                } else if (POT == 'L') {

                    DPhiyDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fy;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);

                // if (DPhiyDrho_old.x != 0 && DPHIDRHO_T[k][i].fy.x == 0) DPHIDRHO_T[k][i].fy.x =1.;
                // if (DPhiyDrho_old.y != 0 && DPHIDRHO_T[k][i].fy.y == 0) DPHIDRHO_T[k][i].fy.y =1.;
                // if (DPhiyDrho_old.z != 0 && DPHIDRHO_T[k][i].fy.z == 0) DPHIDRHO_T[k][i].fy.z =1.;

                // if (DPhiyDrho_old.x != 0. && DPHIDRHO_T[k][i].fy.x == 0. && DPhiyDrho_old.y != 0. && DPHIDRHO_T[k][i].fy.y == 0. && DPhiyDrho_old.z != 0. && DPHIDRHO_T[k][i].fy.z == 0.) {
                //     denom += fabs(DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                // } else {
                //
                //     denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                // }

                denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
                //denom += (DPhiyDrho_old.x*DPhiyDrho_old.x + DPhiyDrho_old.y*DPhiyDrho_old.y + DPhiyDrho_old.z*DPhiyDrho_old.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);

            if (denom != 0.){
                GAMMA[k].y = (SOR)*Phi_old.y/denom;
                GAMMATOT[k].y += GAMMA[k].y;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].y = %.4e\n\n", k, GAMMA[k].y);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];
                    //printf("%.4e \n", GAMMA[k].y);
                    rho_OLD[i].x -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.x; //segmentation fault
                    rho_OLD[i].y -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.y;
                    rho_OLD[i].z -= GAMMA[k].y*DPHIDRHO_T[k][i].fy.z;

                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].y = 0.;
            }

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            //            Z COMPONENT OF THE FORCE
            denom = 0;

            if (POT == 'J') {

                Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z;

            } else if (POT == 'C'){

                Phi_old.z = ShellForce_Cicc(rho_OLD, r_tp1, k).z;

            } else if (POT == 'W'){

                Phi_old.z = Force_WCA(rho_OLD, k).z;

            } else if (POT == 'L'){

                Phi_old.z = Force_LJ(rho_OLD, k).z;

            }

            if(fabs(Phi_old.z)>discr) {

                discr = fabs(Phi_old.z); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.3;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].z = %.4e\n", k, Phi_old.z);

            for (i=0; i<NATOMSPERSPEC[1]; i++) {

                if (POT == 'J') {

                    DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhizDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'W'){

                    //printf("Dsigma_%d_x/Dx_%d = ",k,i);
                    DPhizDrho_old = ConstTens_WCA(rho_OLD, r_tp1, k, i).fz;
                    //printf("(%f %f %f)\n",DPhizDrho_old.x,DPhizDrho_old.y,DPhizDrho_old.z);
                }  else if (POT == 'L'){

                    DPhizDrho_old = ConstTens_LJ(rho_OLD, r_tp1, k, i).fz;

                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);

                // if (DPhizDrho_old.x != 0 && DPHIDRHO_T[k][i].fz.x == 0) DPHIDRHO_T[k][i].fz.x =1.;
                // if (DPhizDrho_old.y != 0 && DPHIDRHO_T[k][i].fz.y == 0) DPHIDRHO_T[k][i].fz.y =1.;
                // if (DPhizDrho_old.z != 0 && DPHIDRHO_T[k][i].fz.z == 0) DPHIDRHO_T[k][i].fz.z =1.;
                //if (k==672)printf("%.8e\n",denom);
                // if (DPhizDrho_old.x != 0. && DPHIDRHO_T[k][i].fz.x == 0. && DPhizDrho_old.y != 0. && DPHIDRHO_T[k][i].fz.y == 0. && DPhizDrho_old.z != 0. && DPHIDRHO_T[k][i].fz.z == 0.) {
                //     denom += fabs(DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                // } else {
                //     denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                // }


                denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
                //denom += (DPhizDrho_old.x*DPhizDrho_old.x + DPhizDrho_old.y*DPhizDrho_old.y + DPhizDrho_old.z*DPhizDrho_old.z);
                // if (k == 672) printf("denom = %.8e\n", denom);
                // if (k == 672) printf("(%d)+= %.8e\n", i,(DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z));
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);
            //if (count == 2 && k == 448) printf("%.4e\n",denom);
            if (denom != 0.){

                GAMMA[k].z = (SOR)*Phi_old.z/denom;
                GAMMATOT[k].z += GAMMA[k].z;

                if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].z = %.4e\n\n", k, GAMMA[k].z);

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    indx_i = INDX[i];

                    rho_OLD[i].x -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.x;
                    rho_OLD[i].y -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.y;
                    rho_OLD[i].z -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.z;

                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho_OLD[i],i);
                }
            }else{
                GAMMA[k].z = 0.;
            }
            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            if (DEBUG_FLAG && _D_CONSTR) printf("it = %d -> Phi[%d] = (%.4e, %.4e, %.4e)\n", count, k, Phi_old.x, Phi_old.y, Phi_old.z);
            //printf("part : %d \n", k);
        } //End loop on constraints
        //if (count == 2) exit(0);
        fprintf(fp_constraints_out, "%d \t %.10e \t %f \n", count, discr, kdiscr);
        testForce = Force_WCA(rho_OLD, 171);
        testdForcex = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fx;
        testdForcey = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fy;
        testdForcez = ConstTens_WCA(rho_OLD, r_tp1, 171, 171).fz;

        printf("nb of iter = %d,\t discr = %e, \t discrk = %.1lf \n", count, discr,kdiscr);
        // printf("indx = %d, shellpos = %.4e %.4e %.4e \n, force = %.4e %.4e %.4e \n",171, rho_OLD[171].x, rho_OLD[171].y, rho_OLD[171].z, testForce.x, testForce.y, testForce.z);
        // printf("dFx = %.4e %.4e %.4e \n dFy = %.4e %.4e %.4e \n dFz = %.4e %.4e %.4e \n", testdForcex.x, testdForcex.y, testdForcex.z, testdForcey.x, testdForcey.y, testdForcey.z, testdForcez.x, testdForcez.y, testdForcez.z);
        // printf("dFxold = %.4e %.4e %.4e \n dFyold = %.4e %.4e %.4e \n dFzold = %.4e %.4e %.4e\n\n", DPHIDRHO_T[171][171].fx.x, DPHIDRHO_T[171][171].fx.y, DPHIDRHO_T[171][171].fx.z, DPHIDRHO_T[171][171].fy.x, DPHIDRHO_T[171][171].fy.y, DPHIDRHO_T[171][171].fy.z, DPHIDRHO_T[171][171].fz.x, DPHIDRHO_T[171][171].fz.y, DPHIDRHO_T[171][171].fz.z);
    } //End while(constraint condition)
    fprintf(fp_constraints_out, "\n");
    fflush(fp_constraints_out);
    fclose(fp_constraints_out);

    SR_ITERS = count;
    SR_DISCR = discr;

//    Uncomment for direct comparison with Conjugate Gradient method
//
//    for (i=0; i<NPART; i++) {
//
//        GAMMA[i] = ShellForce_Jac(rho_OLD, r_tp1, i);
//    }
//
//    SR_DISCR = Variance(GAMMA, NPART);

    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
}

void SteepestDescent(struct point rho[], struct point r[]){

    int i, count = 0;
    double max_Phi = 0, lambda = 2.3e-3, norm = 0., imax_Phi = 0.;
    struct point *Phi_new, *Phi_old, *rho_old, dr, dPhi;

    Phi_new = (struct point *)malloc(NPART*sizeof(struct point));
    Phi_old = (struct point *)malloc(NPART*sizeof(struct point));
    rho_old = (struct point *)malloc(NPART*sizeof(struct point));

    for (i=0; i<NPART; i++) {

        rho_old[i].x = rho[i].x;
        rho_old[i].y = rho[i].y;
        rho_old[i].z = rho[i].z;

        Phi_old[i] = ShellForce_Jac(rho, r, i);

        if (fabs(Phi_old[i].x) > max_Phi) {

            max_Phi = fabs(Phi_old[i].x);
            imax_Phi = i+.1;
        }

        if (fabs(Phi_old[i].y) > max_Phi) {

            max_Phi = fabs(Phi_old[i].y);
            imax_Phi = i+.2;
        }

        if (fabs(Phi_old[i].z) > max_Phi) {

            max_Phi = fabs(Phi_old[i].z);
            imax_Phi = i+.3;
        }
    }

    while(max_Phi > LOW_TOL){

        printf("rho[1] = (%lf, %lf, %lf)\n", rho[1].x, rho[1].y, rho[1].z);

        count++;

        if (count > _MAX_ITER) {

            printf("ShellRelaxation.c -> SteepestDescent() ERROR: Iteration limit exceeded. Convergence not reached!\n");
            exit(EXIT_FAILURE);
        }

        printf("lambda = %lf\n", lambda);

        for (i=0; i<NPART; i++) {

            rho[i].x = rho_old[i].x + lambda*Phi_old[i].x;
            rho[i].y = rho_old[i].y + lambda*Phi_old[i].y;
            rho[i].z = rho_old[i].z + lambda*Phi_old[i].z;
        }

        lambda = 0.;
        norm = 0.;

        for (i=0; i<NPART; i++) {

            Phi_new[i] = ShellForce_Jac(rho, r, i);

            if (fabs(Phi_new[i].x) > max_Phi) {

                max_Phi = fabs(Phi_new[i].x);
                imax_Phi = i+.1;
            }

            if (fabs(Phi_new[i].y) > max_Phi) {

                max_Phi = fabs(Phi_new[i].y);
                imax_Phi = i+.2;
            }

            if (fabs(Phi_new[i].z) > max_Phi) {

                max_Phi = fabs(Phi_new[i].z);
                imax_Phi = i + .3;
            }

            dr.x = rho[i].x - rho_old[i].x;
            dr.y = rho[i].y - rho_old[i].y;
            dr.z = rho[i].z - rho_old[i].z;

            dPhi.x = Phi_new[i].x - Phi_old[i].x;
            dPhi.y = Phi_new[i].y - Phi_old[i].y;
            dPhi.z = Phi_new[i].z - Phi_old[i].z;

            lambda += (dr.x*dPhi.x + dr.y*dPhi.y + dr.z*dPhi.z);
            norm += (dPhi.x*dPhi.x + dPhi.y*dPhi.y + dPhi.z*dPhi.z);

            rho_old[i].x = rho[i].x;
            rho_old[i].y = rho[i].y;
            rho_old[i].z = rho[i].z;

            Phi_old[i].x = Phi_new[i].x;
            Phi_old[i].y = Phi_new[i].y;
            Phi_old[i].z = Phi_new[i].z;
        }

//        lambda /= norm;

        printf("Iteration: %d\tdiscr = (%.25e)\tkdiscr = %.1lf\n", count, max_Phi, imax_Phi);
    }

    free(Phi_new);
    free(Phi_old);
    free(rho_old);
}

void ConjugateGradient(struct point rho[], struct point r[]) {

    int i, line = 0, count = 0;
    double max_Phi = 0.0, denom = 0.0, lambda_down = 0.0, lambda_up = 0.0, lambda1 = 0.0, lambda2 = 0.0, lambda3 = 0.0, beta = 0.0, AA, BB, E0, E1, E2, Emin, sigma = 1.;
    char dir;
    struct point *RHO_OLD, *PHI_OLD, *PHI, *SEARCHDIR;

    RHO_OLD = (struct point *)malloc(NPART*sizeof(struct point));
    PHI_OLD = (struct point *)malloc(NPART*sizeof(struct point));
    PHI = (struct point *)malloc(NPART*sizeof(struct point));
    SEARCHDIR = (struct point *)malloc(NPART*sizeof(struct point));

//    Initialize and check if minimization is required
    for (i=0; i<NPART; i++) {

        RHO_OLD[i].x = rho[i].x;
        RHO_OLD[i].y = rho[i].y;
        RHO_OLD[i].z = rho[i].z;

        PHI_OLD[i] = ShellForce_Jac(rho, r, i);

        if (fabs(PHI_OLD[i].x) > max_Phi) max_Phi = fabs(PHI_OLD[i].x);
        if (fabs(PHI_OLD[i].y) > max_Phi) max_Phi = fabs(PHI_OLD[i].y);
        if (fabs(PHI_OLD[i].z) > max_Phi) max_Phi = fabs(PHI_OLD[i].z);
    }

    if (DEBUG_FLAG && _D_CG) printf("Starting value of max_Phi = %.15e\n", max_Phi);

    SR_DISCR = Variance(PHI_OLD, NPART);

//    Starting minimization process. Looping until ?
    while (SR_DISCR > _SR_CG_LOW_TOL) {

        line++;

        if (line > _MAX_ITER) {

            printf("\nShellRelaxation.c -> ConjugateGradient() ERROR: Iteration limit exceeded. Convergence not reached!\n");
            exit(EXIT_FAILURE);
        }

        beta = 0.0;
        denom = 0.0;
        max_Phi = 0.0;

//        Choosing search direction
        for (i=0; i<NPART; i++) {

            PHI[i] = ShellForce_Jac(rho, r, i);

            if (fabs(PHI[i].x) > max_Phi) max_Phi = fabs(PHI[i].x);
            if (fabs(PHI[i].y) > max_Phi) max_Phi = fabs(PHI[i].y);
            if (fabs(PHI[i].z) > max_Phi) max_Phi = fabs(PHI[i].z);

            denom += (PHI_OLD[i].x*PHI_OLD[i].x + PHI_OLD[i].y*PHI_OLD[i].y + PHI_OLD[i].z*PHI_OLD[i].z);

//            Fletcher & Reeves algorithm
//            beta += (PHI[i].x*PHI[i].x + PHI[i].y*PHI[i].y + PHI[i].z*PHI[i].z);

//            Polak & Ribiere algorithm
            beta += ((PHI[i].x - PHI_OLD[i].x)*PHI[i].x + (PHI[i].y - PHI_OLD[i].y)*PHI[i].y + (PHI[i].z - PHI_OLD[i].z)*PHI[i].z);
        }

//        If first iteration the searching direction is the direction of the force
        (line != 1) ? (beta /= denom) : (beta = 0);

//        Assigning search direction (searchdir)
        for (i=0; i<NPART; i++) {

            SEARCHDIR[i].x = beta*SEARCHDIR[i].x + PHI[i].x;
            SEARCHDIR[i].y = beta*SEARCHDIR[i].y + PHI[i].y;
            SEARCHDIR[i].z = beta*SEARCHDIR[i].z + PHI[i].z;
        }

        if (DEBUG_FLAG && _D_CG) printf("\nNew Search Direction\n");

//        Computing energy before line-minimization process (lambda = 0)
        E0 = EnerPot_Jac(r, rho);
        if (DEBUG_FLAG && _D_CG) printf("E0 = %.4e\n", E0);

        denom = 0.0;
        lambda1 = 0.0;

//        Starting of line minimization process along selected direction
//        Choosing minimization step in the selected direction (lambda)
        for (i=0; i<NPART; i++) {

            denom += K[INDX[i]]*(SEARCHDIR[i].x*SEARCHDIR[i].x + SEARCHDIR[i].y*SEARCHDIR[i].y + SEARCHDIR[i].z*SEARCHDIR[i].z);

            lambda1 += (PHI[i].x*SEARCHDIR[i].x + PHI[i].y*SEARCHDIR[i].y + PHI[i].z*SEARCHDIR[i].z);
        }

//        Saving parameter for following fit (dE/dlambda at lambda = 0)
        BB = -lambda1;
        lambda1 /= denom; //lambda1
        lambda1 *= sigma;

        if (DEBUG_FLAG && _D_CG) printf("denom = %.4e\n", denom);

        for (i=0; i<NPART; i++) {

            rho[i].x = RHO_OLD[i].x + lambda1*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda1*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda1*SEARCHDIR[i].z;
        }

//        Computing energy after minimization process (lambda = lambda1)
        E1 = EnerPot_Jac(r, rho);

//        Fitting parabola to better estimate energy of the minimum
        AA = ((E1-E0) - BB*lambda1)/(lambda1*lambda1);
        if (AA == 0.) printf("\nShellRelaxation.c -> ConjugateGradient() WARNING: NULL AA coefficient: AA = %.10e\n", AA);
//        (AA != 0.) ? (lambda2 = - BB/(2.*AA)) : (lambda2 = lambda1);
        lambda2 = - BB/(2.*AA); //lambda2
        lambda2 *= sigma;

        for (i=0; i<NPART; i++) {

            rho[i].x = RHO_OLD[i].x + lambda2*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda2*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda2*SEARCHDIR[i].z;
        }

//        Computing energy after minimization process and parabola fit (lambda = lambda2)
        E2 = EnerPot_Jac(r, rho);

//        Compute minimum energy (ALERT: CORRECT IF lambda IS ZERO)
        if (E1 < E2){

            Emin = E1;
            lambda_down = lambda1;
            lambda_up = lambda2;

        } else {

            Emin = E2;
            lambda_down = lambda2;
            lambda_up = lambda1;
        }

        for (i=0; i<NPART; i++) {

            rho[i].x = RHO_OLD[i].x + lambda_down*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda_down*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda_down*SEARCHDIR[i].z;
        }

        dir = 'd';

//        Checking minimization success along selected direction
//        If min(E1, E2) < E0 search in selected direction completed. Success and skip to next direction
        if (Emin >= E0) {

//          If min(E1, E2) > E0 and the selected search direction is the first then scale lambda and loop again
            if (line == 1) {

                if (DEBUG_FLAG && _D_CG) printf("iexit == -1 - Reducing lambda and retrying\n");
                if (DEBUG_FLAG && _D_CG) printf("line: %d\tEmin = (%.15e)\tE0 = (%.15e)\n", line, Emin, E0);

                lambda3 = lambda_down*0.1; //lambda3
                lambda3 *= sigma;

                count = 0;

                do {

                    for (i=0; i<NPART; i++) {

                        rho[i].x = RHO_OLD[i].x + lambda3*SEARCHDIR[i].x;
                        rho[i].y = RHO_OLD[i].y + lambda3*SEARCHDIR[i].y;
                        rho[i].z = RHO_OLD[i].z + lambda3*SEARCHDIR[i].z;
                    }

                    Emin = EnerPot_Jac(r, rho);

                    lambda3 /= 2.;

                    count++;

                } while (Emin > E0 || count > _MAX_ITER || lambda3 == 0);

            } else {

//              If min(E1, E2) > E0 the minimum is approaching. Do just another one attempt with lambda scaled.
                lambda3 = lambda2*0.1; //lambda3
                lambda3 *= sigma;

                for (i=0; i<NPART; i++) {

                    rho[i].x = RHO_OLD[i].x + lambda3*SEARCHDIR[i].x;
                    rho[i].y = RHO_OLD[i].y + lambda3*SEARCHDIR[i].y;
                    rho[i].z = RHO_OLD[i].z + lambda3*SEARCHDIR[i].z;
                }

                Emin = EnerPot_Jac(r, rho);

//                If Emin is again higher than E0 than set lambda to be the maximum between lambda1 and lambda2
                if (Emin >= E0) {

                    for (i=0; i<NPART; i++) {

                        rho[i].x = RHO_OLD[i].x + lambda_up*SEARCHDIR[i].x;
                        rho[i].y = RHO_OLD[i].y + lambda_up*SEARCHDIR[i].y;
                        rho[i].z = RHO_OLD[i].z + lambda_up*SEARCHDIR[i].z;
                    }

                    dir = 'u';
                }
            }
        }
//        Terminated search in selected direction

//        Resetting variables after process completion
        for (i=0; i<NPART; i++) {

            RHO_OLD[i].x = rho[i].x;
            RHO_OLD[i].y = rho[i].y;
            RHO_OLD[i].z = rho[i].z;

            PHI_OLD[i].x = PHI[i].x;
            PHI_OLD[i].y = PHI[i].y;
            PHI_OLD[i].z = PHI[i].z;
        }

        for (i=0; i<NPART; i++) {

            PHI[i] = ShellForce_Jac(rho, r, i);

            if (fabs(PHI[i].x) > max_Phi) max_Phi = fabs(PHI[i].x);
            if (fabs(PHI[i].y) > max_Phi) max_Phi = fabs(PHI[i].y);
            if (fabs(PHI[i].z) > max_Phi) max_Phi = fabs(PHI[i].z);
        }

        if (DEBUG_FLAG && _D_CG) printf("max_Phi = %.15e\n", max_Phi);

        SR_DISCR = Variance(PHI, NPART);
        SR_ITERS = line;

        if (SR_DISCR > _UP_TOL) {

            printf("\nShellRelaxation.c -> ConjugateGradient() ERROR: max_Phi = (%.4e). The algorithm has exploded!\n", SR_DISCR);
            exit(EXIT_FAILURE);
        }

        if (line > .1*_MAX_ITER) printf("line = %d\tRMSForce = %.20e\tdir = %c\n", SR_ITERS, SR_DISCR, dir);
    }
}

void MultiConjugateGradient(struct point rho[], struct point r[]) {

    int i, line = 0, count = 0;
    double max_Phi = 0.0, denom = 0.0, lambda_down = 0.0, lambda_up = 0.0, lambda1 = 0.0, lambda2 = 0.0, lambda3 = 0.0, beta = 0.0, AA, BB, E0, E1, E2, Emin, sigma = 1.;
    char dir;
    struct point *RHO_OLD, *PHI_OLD, *PHI, *SEARCHDIR;

    RHO_OLD = (struct point *)malloc(NATOMSPERSPEC[1]*sizeof(struct point));
    PHI_OLD = (struct point *)malloc(NATOMSPERSPEC[1]*sizeof(struct point));
    PHI = (struct point *)malloc(NATOMSPERSPEC[1]*sizeof(struct point));
    SEARCHDIR = (struct point *)malloc(NATOMSPERSPEC[1]*sizeof(struct point));

//    Initialize and check if minimization is required
    for (i=0; i<NATOMSPERSPEC[1]; i++) {

        RHO_OLD[i].x = rho[i].x;
        RHO_OLD[i].y = rho[i].y;
        RHO_OLD[i].z = rho[i].z;

      if (POT == 'W'){

          PHI_OLD[i] = Force_WCA(rho, i);

      }else if (POT == 'L'){

          PHI_OLD[i] = Force_LJ(rho, i);

      }


        if (fabs(PHI_OLD[i].x) > max_Phi) max_Phi = fabs(PHI_OLD[i].x);
        if (fabs(PHI_OLD[i].y) > max_Phi) max_Phi = fabs(PHI_OLD[i].y);
        if (fabs(PHI_OLD[i].z) > max_Phi) max_Phi = fabs(PHI_OLD[i].z);
    }

    if (DEBUG_FLAG && _D_CG) printf("Starting value of max_Phi = %.15e\n", max_Phi);

    SR_DISCR = Variance(PHI_OLD, NATOMSPERSPEC[1]);

//    Starting minimization process. Looping until ?
    while (SR_DISCR > _SR_CG_LOW_TOL) {
        printf("%.4e \n", SR_DISCR);
        line++;

        if (line > _MAX_ITER) {

            printf("\nShellRelaxation.c -> ConjugateGradient() ERROR: Iteration limit exceeded. Convergence not reached!\n");
            exit(EXIT_FAILURE);
        }

        beta = 0.0;
        denom = 0.0;
        max_Phi = 0.0;

//        Choosing search direction
        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            if (POT == 'W'){

                PHI[i] = Force_WCA(rho, i);

            }else if (POT == 'L'){

                PHI[i] = Force_LJ(rho, i);

            }

            if (fabs(PHI[i].x) > max_Phi) max_Phi = fabs(PHI[i].x);
            if (fabs(PHI[i].y) > max_Phi) max_Phi = fabs(PHI[i].y);
            if (fabs(PHI[i].z) > max_Phi) max_Phi = fabs(PHI[i].z);

            denom += (PHI_OLD[i].x*PHI_OLD[i].x + PHI_OLD[i].y*PHI_OLD[i].y + PHI_OLD[i].z*PHI_OLD[i].z);

//            Fletcher & Reeves algorithm
//            beta += (PHI[i].x*PHI[i].x + PHI[i].y*PHI[i].y + PHI[i].z*PHI[i].z);

//            Polak & Ribiere algorithm
            beta += ((PHI[i].x - PHI_OLD[i].x)*PHI[i].x + (PHI[i].y - PHI_OLD[i].y)*PHI[i].y + (PHI[i].z - PHI_OLD[i].z)*PHI[i].z);
        }

//        If first iteration the searching direction is the direction of the force
        (line != 1) ? (beta /= denom) : (beta = 0);

//        Assigning search direction (searchdir)
        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            SEARCHDIR[i].x = beta*SEARCHDIR[i].x + PHI[i].x;
            SEARCHDIR[i].y = beta*SEARCHDIR[i].y + PHI[i].y;
            SEARCHDIR[i].z = beta*SEARCHDIR[i].z + PHI[i].z;
        }

        if (DEBUG_FLAG && _D_CG) printf("\nNew Search Direction\n");

//        Computing energy before line-minimization process (lambda = 0)
        if (POT == 'W'){

            E0 = EnerPot_WCA(rho);

        }else if (POT == 'L'){

            E0 = EnerPot_LJ(rho);

        }
        if (DEBUG_FLAG && _D_CG) printf("E0 = %.4e\n", E0);

        denom = 0.0;
        lambda1 = 0.0;

//        Starting of line minimization process along selected direction
//        Choosing minimization step in the selected direction (lambda)
        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            denom += (SEARCHDIR[i].x*SEARCHDIR[i].x + SEARCHDIR[i].y*SEARCHDIR[i].y + SEARCHDIR[i].z*SEARCHDIR[i].z);//%*K[INDX[i]]?;

            lambda1 += (PHI[i].x*SEARCHDIR[i].x + PHI[i].y*SEARCHDIR[i].y + PHI[i].z*SEARCHDIR[i].z);
        }

//        Saving parameter for following fit (dE/dlambda at lambda = 0)
        BB = -lambda1;
        lambda1 /= denom; //lambda1
        lambda1 *= sigma;

        if (DEBUG_FLAG && _D_CG) printf("denom = %.4e\n", denom);

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            rho[i].x = RHO_OLD[i].x + lambda1*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda1*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda1*SEARCHDIR[i].z;
            Rem_Point_From_Cell(i);
            Add_Point_To_Cell(rho[i],i);
        }

//        Computing energy after minimization process (lambda = lambda1)
        if (POT == 'W'){

            E1 = EnerPot_WCA(rho);

        }else if (POT == 'L'){

            E1 = EnerPot_LJ(rho);

        }

//        Fitting parabola to better estimate energy of the minimum
        AA = ((E1-E0) - BB*lambda1)/(lambda1*lambda1);
        if (AA == 0.) printf("\nShellRelaxation.c -> ConjugateGradient() WARNING: NULL AA coefficient: AA = %.10e\n", AA);
//        (AA != 0.) ? (lambda2 = - BB/(2.*AA)) : (lambda2 = lambda1);
        lambda2 = - BB/(2.*AA); //lambda2
        lambda2 *= sigma;

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            rho[i].x = RHO_OLD[i].x + lambda2*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda2*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda2*SEARCHDIR[i].z;
            Rem_Point_From_Cell(i);
            Add_Point_To_Cell(rho[i],i);
        }

//        Computing energy after minimization process and parabola fit (lambda = lambda2)
        if (POT == 'W'){

            E2 = EnerPot_WCA(rho);

        }else if (POT == 'L'){

            E2 = EnerPot_LJ(rho);

        }

//        Compute minimum energy (ALERT: CORRECT IF lambda IS ZERO)
        if (E1 < E2){

            Emin = E1;
            lambda_down = lambda1;
            lambda_up = lambda2;

        } else {

            Emin = E2;
            lambda_down = lambda2;
            lambda_up = lambda1;
        }

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            rho[i].x = RHO_OLD[i].x + lambda_down*SEARCHDIR[i].x;
            rho[i].y = RHO_OLD[i].y + lambda_down*SEARCHDIR[i].y;
            rho[i].z = RHO_OLD[i].z + lambda_down*SEARCHDIR[i].z;
            Rem_Point_From_Cell(i);
            Add_Point_To_Cell(rho[i],i);
        }

        dir = 'd';

//        Checking minimization success along selected direction
//        If min(E1, E2) < E0 search in selected direction completed. Success and skip to next direction
        if (Emin >= E0) {

//          If min(E1, E2) > E0 and the selected search direction is the first then scale lambda and loop again
            if (line == 1) {

                if (DEBUG_FLAG && _D_CG) printf("iexit == -1 - Reducing lambda and retrying\n");
                if (DEBUG_FLAG && _D_CG) printf("line: %d\tEmin = (%.15e)\tE0 = (%.15e)\n", line, Emin, E0);

                lambda3 = lambda_down*0.1; //lambda3
                lambda3 *= sigma;

                count = 0;

                do {

                    for (i=0; i<NATOMSPERSPEC[1]; i++) {

                        rho[i].x = RHO_OLD[i].x + lambda3*SEARCHDIR[i].x;
                        rho[i].y = RHO_OLD[i].y + lambda3*SEARCHDIR[i].y;
                        rho[i].z = RHO_OLD[i].z + lambda3*SEARCHDIR[i].z;
                        Rem_Point_From_Cell(i);
                        Add_Point_To_Cell(rho[i],i);
                    }

                    if (POT == 'W'){

                        Emin = EnerPot_WCA(rho);

                    }else if (POT == 'L'){

                        Emin = EnerPot_LJ(rho);

                    }

                    lambda3 /= 2.;

                    count++;

                } while (Emin > E0 || count > _MAX_ITER || lambda3 == 0);

            } else {

//              If min(E1, E2) > E0 the minimum is approaching. Do just another one attempt with lambda scaled.
                lambda3 = lambda2*0.1; //lambda3
                lambda3 *= sigma;

                for (i=0; i<NATOMSPERSPEC[1]; i++) {

                    rho[i].x = RHO_OLD[i].x + lambda3*SEARCHDIR[i].x;
                    rho[i].y = RHO_OLD[i].y + lambda3*SEARCHDIR[i].y;
                    rho[i].z = RHO_OLD[i].z + lambda3*SEARCHDIR[i].z;
                    Rem_Point_From_Cell(i);
                    Add_Point_To_Cell(rho[i],i);
                }

                if (POT == 'W'){

                    Emin = EnerPot_WCA(rho);

                }else if (POT == 'L'){

                    Emin = EnerPot_LJ(rho);

                }

//                If Emin is again higher than E0 than set lambda to be the maximum between lambda1 and lambda2
                if (Emin >= E0) {

                    for (i=0; i<NATOMSPERSPEC[1]; i++) {

                        rho[i].x = RHO_OLD[i].x + lambda_up*SEARCHDIR[i].x;
                        rho[i].y = RHO_OLD[i].y + lambda_up*SEARCHDIR[i].y;
                        rho[i].z = RHO_OLD[i].z + lambda_up*SEARCHDIR[i].z;
                        Rem_Point_From_Cell(i);
                        Add_Point_To_Cell(rho[i],i);
                    }

                    dir = 'u';
                }
            }
        }
//        Terminated search in selected direction

//        Resetting variables after process completion
        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            RHO_OLD[i].x = rho[i].x;
            RHO_OLD[i].y = rho[i].y;
            RHO_OLD[i].z = rho[i].z;

            PHI_OLD[i].x = PHI[i].x;
            PHI_OLD[i].y = PHI[i].y;
            PHI_OLD[i].z = PHI[i].z;
        }

        for (i=0; i<NATOMSPERSPEC[1]; i++) {

            if (POT == 'W'){

                PHI[i] = Force_WCA(rho, i);

            }else if (POT == 'L'){

                PHI[i] = Force_LJ(rho, i);

            }

            if (fabs(PHI[i].x) > max_Phi) max_Phi = fabs(PHI[i].x);
            if (fabs(PHI[i].y) > max_Phi) max_Phi = fabs(PHI[i].y);
            if (fabs(PHI[i].z) > max_Phi) max_Phi = fabs(PHI[i].z);
            //printf("F_%d = %.4e %.4e %.4e\n", 0, PHI[0].x, PHI[0].y, PHI[0].z);
        }

        if (DEBUG_FLAG && _D_CG) printf("max_Phi = %.15e\n", max_Phi);

        SR_DISCR = Variance(PHI, NATOMSPERSPEC[1]);
        SR_ITERS = line;

        if (SR_DISCR > _UP_TOL) {

            printf("\nShellRelaxation.c -> ConjugateGradient() ERROR: max_Phi = (%.4e). The algorithm has exploded!\n", SR_DISCR);
            exit(EXIT_FAILURE);
        }

        if (line > .1*_MAX_ITER) printf("line = %d\tRMSForce = %.20e\tdir = %c\n", SR_ITERS, SR_DISCR, dir);
    }
}
