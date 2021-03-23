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

    //    struct point s;

    for (k=0; k<NPART; k++) {

        if (POT == 'J') {

            Phi_old = ShellForce_Jac(rho_OLD, r_tp1, k);

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

                DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)

            }else if (POT == 'C') {

                DPHIDRHO_T[k][i] = ConstTens_Cicc(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)
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

            } else if (POT == 'C'){

                Phi_old.x = ShellForce_Cicc(rho_OLD, r_tp1, k).x;
            }

            if(fabs(Phi_old.x)>discr) {

                discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
                kdiscr = k + 0.1;
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);

            for (i=0; i<NPART; i++) {

                if (POT == 'J') {

                    DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD

                } else if (POT == 'C'){

                    DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD
                }

                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
            }

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);

            GAMMA[k].x = Phi_old.x/denom;

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
    } //End while(constraint condition)

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

void BSHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], struct point vrho_t[], struct point vrho_OLD[], int timestep, int ccount){

    int k, i, indx_i, count = 0;
    double discr = 0, kdiscr = -1.;
    double denom;
    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old;
    struct point DPhixDvrho_old, DPhiyDvrho_old;

    //    struct point s;

    for (k=0; k<NPART; k++) {

        if (POT == 'J') {

            Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x + 0.5*CHI[INDX[k]]*B0*vrho_OLD[k].y;//predicted?
            Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y - 0.5*CHI[INDX[k]]*B0*vrho_OLD[k].x;
            Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z;
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
        DPHIDVRHO_T[k][k].fx.y = 0.5*CHI[INDX[k]]*B0;

        DPHIDVRHO_T[k][k].fy.x = -0.5*CHI[INDX[k]]*B0;
    }
    printf("discr = %e\n", discr);

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

        for (k=0; k<NPART; k++){

            SHELLACC_TM1[k].x = 0;
            SHELLACC_TM1[k].y = 0;
            SHELLACC_TM1[k].z = 0;
        }

        for (k=0; k<NPART; k++) { //Looping on all constraints

            //            X COMPONENT OF THE FORCE
            //denom = 0;

            if (POT == 'J') {

                Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x + 0.5*CHI[INDX[k]]*B0*vrho_OLD[k].y;
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

                DPhixDvrho_old.y = 0.5*CHI[INDX[k]]*B0;


            } else if (POT == 'C'){

                DPhixDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, k).fx; // TO BE COMPUTED FOR r_OLD
            }

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);

            // if (B0 == 0){
            //     denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
            // }

                //denom += DPhixDrho_old.y*DPHIDVRHO_T[k][i].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][i].fx.y;
            denom = 0.5*DT*(DPhixDrho_old.y*DPHIDVRHO_T[k][k].fx.y*DT + DPhixDvrho_old.y*DPHIDVRHO_T[k][k].fx.y*3);
            // printf("denom.x (part %d) = %.4e \t",denom, k);

            //}

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);


            GAMMA[k].x = Phi_old.x/denom;
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

                Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y - 0.5*CHI[INDX[k]]*B0*vrho_OLD[k].x;
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

                DPhiyDvrho_old.x = -0.5*CHI[INDX[k]]*B0;

            } else if (POT == 'C'){

                DPhiyDrho_old = ConstTens_Cicc(rho_OLD, r_tp1, k, k).fy; // TO BE COMPUTED FOR r_OLD
            }

            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);

            //denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);

            denom = 0.5*DT*(DPhiyDrho_old.x*DPHIDVRHO_T[k][k].fy.x*DT + DPhiyDvrho_old.x*DPHIDVRHO_T[k][k].fy.x*3);
            // printf("denom.y (part %d) = %.4e \t",denom, k);


            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);

            GAMMA[k].y = Phi_old.y/denom;
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
                denom += 0.5*DT*(DPhizDrho_old.x*DPHIDVRHO_T[k][i].fz.x*DT + DPhizDrho_old.y*DPHIDVRHO_T[k][i].fz.y*DT + DPhizDrho_old.z*DPHIDVRHO_T[k][i].fz.z*DT);
            }
            // printf("denom.z (part %d) = %.4e \n",denom, k);

            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);

            GAMMA[k].z = Phi_old.z/denom;

            if (DEBUG_FLAG && _D_SHAKE) printf("GAMMA[%d].z = %.4e\n\n", k, GAMMA[k].z);

            // for (i=0; i<NPART; i++) {
            //
            //     indx_i = INDX[i];
            //
            //     rho_OLD[i].x -= GAMMA[k].x*DPHIDRHO_T[k][i].fx.x;
            //     rho_OLD[i].y -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.y;
            //     rho_OLD[i].z -= GAMMA[k].z*DPHIDRHO_T[k][i].fz.z;
            // }
            rho_OLD[k].x -= 0.5*DT*DT*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            rho_OLD[k].y -= 0.5*DT*DT*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            rho_OLD[k].z -= 0;

            vrho_OLD[k].x -= (1.5*DT*GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x);
            vrho_OLD[k].y -= (1.5*DT*GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y);
            vrho_OLD[k].z -= 0;

            // SHELLACC_TM1[k].x += GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x;
            // SHELLACC_TM1[k].y += GAMMA[k].x*DPHIDVRHO_T[k][k].fx.y;
            // SHELLACC_TM1[k].z += 0;

            for (i=0; i<NPART; i++) {

                rho_OLD[i].x -= (0.5*DT*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
                rho_OLD[i].y -= (0.5*DT*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
                rho_OLD[i].z -= (0.5*DT*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);

                vrho_OLD[i].x -= (1.5*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.x);
                vrho_OLD[i].y -= (1.5*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.y);
                vrho_OLD[i].z -= (1.5*DT*GAMMA[k].z*DPHIDVRHO_T[k][i].fz.z);

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
            }

            // printf("Shell acc (part %d , iter %d) = (%.4e,%.4e,%.4e) \t",k,count,SHELLACC_TM1[k].x,SHELLACC_TM1[k].y,SHELLACC_TM1[k].z);
            // printf("Phi (part %d , iter %d) = (%.4e,%.4e,%.4e) \n",k,count, Phi_old.x, Phi_old.y, Phi_old.z);
            // rho_OLD[k].x -= (GAMMA[k].x*DPHIDVRHO_T[k][k].fx.x+GAMMA[k].y*DPHIDVRHO_T[k][k].fy.x)*0.5*DT*DT;

            if (DEBUG_FLAG && _D_SHAKE)  {

                for (i=0; i<NPART; i++) {

                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
                }
            }

            if (DEBUG_FLAG && _D_CONSTR) printf("it = %d -> Phi[%d] = (%.4e, %.4e, %.4e)\n", count, k, Phi_old.x, Phi_old.y, Phi_old.z);
        } //End loop on constraints
        //printf("nb of iter = %d,\t discr = %e, \t discrk = %.1lf \n", count, discr,kdiscr);

    } //End while(constraint condition)
    // for (i=0; i<NPART; i++) {
    //
    //     SHELLACC_TM1[i].x = GAMMA[i].x*DPHIDVRHO_T[i][i].fx.x+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.x+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.x;
    //     SHELLACC_TM1[i].y = GAMMA[i].x*DPHIDVRHO_T[i][i].fx.y+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.y+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.y;
    //     SHELLACC_TM1[i].z = GAMMA[i].x*DPHIDVRHO_T[i][i].fx.z+GAMMA[i].y*DPHIDVRHO_T[i][i].fy.z+GAMMA[i].z*DPHIDVRHO_T[i][i].fz.z;
    // }

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
