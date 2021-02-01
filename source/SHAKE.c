//
//  ShaPMoD.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "SHAKE.h"

void ML_SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount){
    
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
    
    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
}

//void SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[]){
//
//    int k, i, indx_i, count = 0;
//    double discr = 0, gamma = 0;
//    double denom;
//    struct point Phi_old, DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old;
//
//    for (k=0; k<NPART; k++) {
//
//        Phi_old = ShellForce_Jac(rho_OLD, r_tp1, k);
//
//        if (fabs(Phi_old.x) > discr) discr = fabs(Phi_old.x);
//        if (fabs(Phi_old.y) > discr) discr = fabs(Phi_old.y);
//        if (fabs(Phi_old.z) > discr) discr = fabs(Phi_old.z);
//
//        for (i=0; i<NPART; i++) {
//
//            DPHIDRHO_T[k][i] = ConstTens_Jac(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t)
//
//            if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) {
//                printf("DPHIDRHO_T[%d][%d] =\n", k, i);
//                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", DPHIDRHO_T[k][i].fx.x, DPHIDRHO_T[k][i].fx.y, DPHIDRHO_T[k][i].fx.z, DPHIDRHO_T[k][i].fy.x, DPHIDRHO_T[k][i].fy.y, DPHIDRHO_T[k][i].fy.z, DPHIDRHO_T[k][i].fz.x, DPHIDRHO_T[k][i].fz.y, DPHIDRHO_T[k][i].fz.z);
//            }
//        }
//    }
//
//    while (discr > LOW_TOL) { //Verifying the constraint condition
//
//        if (VERBOSE_FLAG && _V_SHAKE){
//
//            printf("Iteration: %d\tdiscr = (%.4e)\n", count, discr);
//        }
//
//        Write_SHAKE_output(1, count, discr);
//
//        count++;
//
//        if (count>_MAX_ITER) {
//
//            printf("\nSHAKE.c -> SHAKE ERROR: Iteration time limit exceeded. Convergence not reached!\n");
//            exit(EXIT_FAILURE);
//
//        }else if (discr>_UP_TOL){
//
//            printf("\nSHAKE.c -> SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
//            exit(EXIT_FAILURE);
//        }
//
//        discr = 0;
//
//        for (k=0; k<NPART; k++) { //Looping on all constraints
//
//            //            X COMPONENT OF THE FORCE
//            denom = 0;
//
//            if((fabs((Phi_old.x = ShellForce_Jac(rho_OLD, r_tp1, k).x)))>discr) discr = fabs(Phi_old.x); // TO BE COMPUTED FOR r_OLD
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("Phix_old[%d].x = %.4e\n", k, Phi_old.x);
//
//            for (i=0; i<NPART; i++) {
//
//                DPhixDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD
//                
//                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhixDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhixDrho_old.x, DPhixDrho_old.y, DPhixDrho_old.z);
//
//                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z)/MU[INDX[i]];
//            }
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fx = %.4e\n", k, k, denom);
//
//            gamma = Phi_old.x/denom;
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("gamma[%d] = %.4e\n\n", k, gamma);
//
//            for (i=0; i<NPART; i++) {
//
//                indx_i = INDX[i];
//
//                rho_OLD[i].x -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fx.x;
//                rho_OLD[i].y -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fx.y;
//                rho_OLD[i].z -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fx.z;
//            }
//
//            if (DEBUG_FLAG && _D_SHAKE)  {
//
//                for (i=0; i<NPART; i++) {
//
//                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
//                }
//            }
//
//            //            Y COMPONENT OF THE FORCE
//            denom = 0;
//            
//            if((fabs((Phi_old.y = ShellForce_Jac(rho_OLD, r_tp1, k).y)))>discr) discr = fabs(Phi_old.y); // TO BE COMPUTED FOR r_OLD
//            
//            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].y = %.4e\n", k, Phi_old.y);
//
//            for (i=0; i<NPART; i++) {
//
//                DPhiyDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD
//
//                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhiyDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhiyDrho_old.x, DPhiyDrho_old.y, DPhiyDrho_old.z);
//
//                denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z)/MU[INDX[i]];
//            }
//            
//            if (DEBUG_FLAG && _D_SHAKE) printf("DPhiyDrho_old[%d] dot DPHIDRHO_T[%d].fy = %.4e\n", k, k, denom);
//
//            gamma = Phi_old.y/denom;
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("gamma[%d] = %.4e\n\n", k, gamma);
//
//            for (i=0; i<NPART; i++) {
//
//                indx_i = INDX[i];
//
//                rho_OLD[i].x -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fy.x;
//                rho_OLD[i].y -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fy.y;
//                rho_OLD[i].z -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fy.z;
//            }
//
//            if (DEBUG_FLAG && _D_SHAKE)  {
//
//                for (i=0; i<NPART; i++) {
//
//                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
//                }
//            }
//
//            //            Z COMPONENT OF THE FORCE
//            denom = 0;
//
//            if((fabs((Phi_old.z = ShellForce_Jac(rho_OLD, r_tp1, k).z)))>discr) discr = fabs(Phi_old.z); // TO BE COMPUTED FOR r_OLD
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("sigma_old[%d].z = %.4e\n", k, Phi_old.z);
//
//            for (i=0; i<NPART; i++) {
//
//                DPhizDrho_old = ConstTens_Jac(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD
//
//                if (DEBUG_FLAG && _D_SHAKE && _D_TENSOR) printf("DPhizDrho_old[%d][%d] = (%.4e, %.4e, %.4e)\n", k, i, DPhizDrho_old.x, DPhizDrho_old.y, DPhizDrho_old.z);
//
//                denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z)/MU[INDX[i]];
//            }
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("DPhizDrho_old[%d] dot DPHIDRHO_T[%d].fz = %.4e\n", k, k, denom);
//
//            gamma = Phi_old.z/denom;
//
//            if (DEBUG_FLAG && _D_SHAKE) printf("gamma[%d] = %.4e\n\n", k, gamma);
//
//            for (i=0; i<NPART; i++) {
//
//                indx_i = INDX[i];
//
//                rho_OLD[i].x -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fz.x;
//                rho_OLD[i].y -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fz.y;
//                rho_OLD[i].z -= gamma/MU[indx_i]*DPHIDRHO_T[k][i].fz.z;
//            }
//
//            if (DEBUG_FLAG && _D_SHAKE)  {
//
//                for (i=0; i<NPART; i++) {
//
//                    printf("rho_NEW[%d] = (%.4e, %.4e, %.4e)\n", i, rho_OLD[i].x, rho_OLD[i].y, rho_OLD[i].z);
//                }
//            }
//        } //End loop on constraints
//    } //End while(constraint condition)
//
//    SR_ITERS = count;
//    SR_DISCR = discr;
//
//    if (VERBOSE_FLAG && _V_SHAKE) printf("Convergence of SHAKE reached after %d iteration(s)\ndiscr = (%.4e)\n\n", count, discr);
//}

//void ML_SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[]){
//
//    int k, i, count = 0;
//    double discr = {0}, gamma = {0};
//    double sigma_old, denom;
//    struct point DPhixDrho_old, DPhiyDrho_old, DPhizDrho_old;
//
//    for (k=0; k<_NPART; k++) {
//
//        for (i=0; i<_NPART; i++) {
//
//            DPHIDRHO_T[k][i] = ConstTens(rho_t, r_t, k, i); // TO BE COMPUTED FOR r(t) = SHELLPOS_T
//        }
//    }
//
//    do {
//
//        discr = 0;
//
//        for (k=0; k<_NPART; k++) {
//
//            //            X COMPONENT OF THE FORCE
//            denom = 0;
//
//            if((fabs((sigma_old = ShellForce(rho_OLD, r_tp1, k).x)))>discr) discr = fabs(sigma_old);; // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//            for (i=0; i<_NPART; i++) {
//
//                DPhixDrho_old = ConstTens(rho_OLD, r_tp1, k, i).fx; // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//                denom += (DPhixDrho_old.x*DPHIDRHO_T[k][i].fx.x + DPhixDrho_old.y*DPHIDRHO_T[k][i].fx.y + DPhixDrho_old.z*DPHIDRHO_T[k][i].fx.z);
//            }
//
//            gamma = sigma_old/denom;
//
//            for (i=0; i<_NPART; i++) {
//
//                rho_OLD[i].x -= gamma*DPHIDRHO_T[k][i].fx.x;
//                rho_OLD[i].y -= gamma*DPHIDRHO_T[k][i].fx.y;
//                rho_OLD[i].z -= gamma*DPHIDRHO_T[k][i].fx.z;
//            }
//
//            //            Y COMPONENT OF THE FORCE
//            denom = 0;
//
//            if((fabs((sigma_old = ShellForce(rho_OLD, r_tp1, k).y)))>discr) discr = fabs(sigma_old); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//            for (i=0; i<_NPART; i++) {
//
//                DPhiyDrho_old = ConstTens(rho_OLD, r_tp1, k, i).fy; // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//                denom += (DPhiyDrho_old.x*DPHIDRHO_T[k][i].fy.x + DPhiyDrho_old.y*DPHIDRHO_T[k][i].fy.y + DPhiyDrho_old.z*DPHIDRHO_T[k][i].fy.z);
//            }
//
//            gamma = sigma_old/denom;
//
//            for (i=0; i<_NPART; i++) {
//
//                rho_OLD[i].x -= gamma*DPHIDRHO_T[k][i].fy.x;
//                rho_OLD[i].y -= gamma*DPHIDRHO_T[k][i].fy.y;
//                rho_OLD[i].z -= gamma*DPHIDRHO_T[k][i].fy.z;
//            }
//
//            //            Z COMPONENT OF THE FORCE
//            denom = 0;
//
//            if((fabs((sigma_old = ShellForce(rho_OLD, r_tp1, k).z)))>discr) discr = fabs(sigma_old); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//            for (i=0; i<_NPART; i++) {
//
//                DPhizDrho_old = ConstTens(rho_OLD, r_tp1, k, i).fz; // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//                denom += (DPhizDrho_old.x*DPHIDRHO_T[k][i].fz.x + DPhizDrho_old.y*DPHIDRHO_T[k][i].fz.y + DPhizDrho_old.z*DPHIDRHO_T[k][i].fz.z);
//            }
//
//            gamma = sigma_old/denom;
//
//            for (i=0; i<_NPART; i++) {
//
//                rho_OLD[i].x -= gamma*DPHIDRHO_T[k][i].fz.x;
//                rho_OLD[i].y -= gamma*DPHIDRHO_T[k][i].fz.y;
//                rho_OLD[i].z -= gamma*DPHIDRHO_T[k][i].fz.z;
//            }
//        }
//
//        count++;
//
//        if (count>_MAX_ITER) {
//
//            printf("SHAKE.c -> ML_SHAKE ERROR: Iteration time limit exceeded. Convergence not reached!\n");
//            exit(EXIT_FAILURE);
//
//        }else if (discr>_UP_TOL){
//
//            printf("SHAKE.c -> ML_SHAKE ERROR: discr = (%.4e). The algorithm has exploded!\n", discr);
//            exit(EXIT_FAILURE);
//
//        }else if (_DEBUG && _D_SHAKE){
//
//            printf("Iteration: %d\tdiscr = (%.4e)\n", count, discr);
//        }
//
//    } while (discr>_LOW_TOL);
//
//    printf("Convergence of SHAKE reached after %d iteration(s)\n", count);
//    printf("discr = (%.4e)\n\n", discr);
//}
//
//void VEC_SHAKE(struct point r_t[], struct point r_OLD[]){
//
//    int k, i, count = 0;
//    double DT2overMU = (_DT*_DT/_MU);
//    struct point maxsigma = {0}, gamma = {0};
//    struct point sigma_old, tensdottens;
//    struct tensor W_t[_NPART][_NPART] = {0}, W_old = {0};
//
//    for (k=0; k<_NPART; k++) {
//
//        gamma.x = 0;
//        gamma.y = 0;
//        gamma.z = 0;
//
//        for (i=0; i<_NPART; i++) {
//
//            W_t[k][i] = ConstTens(r_t, k, i); // TO BE COMPUTED FOR r(t) = SHELLPOS_T
//
//            if (_DEBUG && _D_SHAKE && _D_TENSOR) {
//                printf("W_t[%d][%d] =\n", k, i);
//                printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", W_t[k][i].fx.x, W_t[k][i].fx.y, W_t[k][i].fx.z, W_t[k][i].fy.x, W_t[k][i].fy.y, W_t[k][i].fy.z, W_t[k][i].fz.x, W_t[k][i].fz.y, W_t[k][i].fz.z);
//            }
//        }
//    }
//
//    do {
//
//        maxsigma.x = maxsigma.y = maxsigma.z = 0;
//
//        for (k=0; k<_NPART; k++) {
//
//            tensdottens.x = 0;
//            tensdottens.y = 0;
//            tensdottens.z = 0;
//
//            sigma_old = ShellForce(r_OLD, PARTPOS_TP1, k); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//            if (maxsigma.x<fabs(sigma_old.x)) maxsigma.x = fabs(sigma_old.x);
//            if (maxsigma.y<fabs(sigma_old.y)) maxsigma.y = fabs(sigma_old.y);
//            if (maxsigma.z<fabs(sigma_old.z)) maxsigma.z = fabs(sigma_old.z);
//
//            if (_DEBUG && _D_SHAKE) printf("sigma_old[%d] = (%.4e, %.4e, %.4e)\n", k, sigma_old.x, sigma_old.y, sigma_old.z);
//
//            for (i=0; i<_NPART; i++) {
//
//                W_old = ConstTens(r_OLD, k, i); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//                if (_DEBUG && _D_SHAKE && _D_TENSOR)  {
//                    printf("W_old[%d][%d] =\n", k, i);
//                    printf("%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n%.4e\t%.4e\t%.4e\n\n", W_old.fx.x, W_old.fx.y, W_old.fx.z, W_old.fy.x, W_old.fy.y, W_old.fy.z, W_old.fz.x, W_old.fz.y, W_old.fz.z);
//                }
//
//                tensdottens.x += W_old.fx.x*W_t[k][i].fx.x + W_old.fx.y*W_t[k][i].fx.y + W_old.fx.z*W_t[k][i].fx.z;
//                tensdottens.y += W_old.fy.x*W_t[k][i].fy.x + W_old.fy.y*W_t[k][i].fy.y + W_old.fy.z*W_t[k][i].fy.z;
//                tensdottens.z += W_old.fz.x*W_t[k][i].fz.x + W_old.fz.y*W_t[k][i].fz.y + W_old.fz.z*W_t[k][i].fz.z;
//            }
//
//            if (_DEBUG && _D_SHAKE) printf("W_old[%d] dot W_t[%d] = (%.4e, %.4e, %.4e)\n", k, k, tensdottens.x, tensdottens.y, tensdottens.z);
//
//            gamma.x = sigma_old.x/DT2overMU/tensdottens.x;
//            gamma.y = sigma_old.y/DT2overMU/tensdottens.y;
//            gamma.z = sigma_old.z/DT2overMU/tensdottens.z;
//
//            if (_DEBUG && _D_SHAKE) printf("gamma[%d] = (%.4e, %.4e, %.4e)\n\n", k, gamma.x, gamma.y, gamma.z);
//
//            for (i=0; i<_NPART; i++) {
//
//                r_OLD[i].x -= DT2overMU*(gamma.x*W_t[k][i].fx.x + gamma.y*W_t[k][i].fy.x + gamma.z*W_t[k][i].fz.x);
//                r_OLD[i].y -= DT2overMU*(gamma.x*W_t[k][i].fx.y + gamma.y*W_t[k][i].fy.y + gamma.z*W_t[k][i].fz.y);
//                r_OLD[i].z -= DT2overMU*(gamma.x*W_t[k][i].fx.z + gamma.y*W_t[k][i].fy.z + gamma.z*W_t[k][i].fz.z);
//            }
//        }
//
//        if (_DEBUG && _D_SHAKE)  {
//
//            for (i=0; i<_NPART; i++) {
//
//                printf("r_OLD[%d] = (%.4e, %.4e, %.4e)\n", i, r_OLD[i].x, r_OLD[i].y, r_OLD[i].z);
//            }
//        }
//
//        count++;
//
//        if (count>_MAX_ITER) {
//
//            printf("SHAKE.c -> ERROR: Iteration time limit exceeded. Convergence not reached!!\n");
//            exit(EXIT_FAILURE);
//
//        }else{
//
//            if (_DEBUG && _D_SHAKE) printf("Iteration: %d\tmaxsigma = (%.4e, %.4e, %.4e)\n", count, maxsigma.x, maxsigma.y, maxsigma.z);
//        }
//
//    } while (maxsigma.x>_TOL || maxsigma.y>_TOL || maxsigma.z>_TOL);
//
//    printf("Convergence of SHAKE reached after %d iteration(s)\n", count);
//    printf("maxsigma = (%.4e, %.4e, %.4e)\n\n", maxsigma.x, maxsigma.y, maxsigma.z);
//}
//
//void ML_VEC_SHAKE(struct point r_t[], struct point r_OLD[]){
//
//    int k, i, count = 0;
//    double DT2 = _DT*_DT;
//    struct point maxsigma = {0}, gamma = {0};
//    struct point sigma_old, tensdottens;
//    struct tensor W_t[_NPART][_NPART] = {0}, W_old = {0};
//
//    for (k=0; k<_NPART; k++) {
//
//        gamma.x = 0;
//        gamma.y = 0;
//        gamma.z = 0;
//
//        for (i=0; i<_NPART; i++) {
//
//            W_t[k][i] = ConstTens(r_t, k, i); // TO BE COMPUTED FOR r(t) = SHELLPOS_T
//        }
//    }
//
//    do {
//
//        maxsigma.x = maxsigma.y = maxsigma.z = 0;
//
//        for (k=0; k<_NPART; k++) {
//
//            tensdottens.x = 0;
//            tensdottens.y = 0;
//            tensdottens.z = 0;
//
//            sigma_old = ShellForce(r_OLD, PARTPOS_TP1, k); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//            if (maxsigma.x<fabs(sigma_old.x)) maxsigma.x = fabs(sigma_old.x);
//            if (maxsigma.y<fabs(sigma_old.y)) maxsigma.y = fabs(sigma_old.y);
//            if (maxsigma.z<fabs(sigma_old.z)) maxsigma.z = fabs(sigma_old.z);
//
//            for (i=0; i<_NPART; i++) {
//
//                W_old = ConstTens(r_OLD, k, i); // TO BE COMPUTED FOR r_OLD = SHELLPOS_Tt
//
//                tensdottens.x += W_old.fx.x*W_t[k][i].fx.x + W_old.fx.y*W_t[k][i].fx.y + W_old.fx.z*W_t[k][i].fx.z;
//                tensdottens.y += W_old.fy.x*W_t[k][i].fy.x + W_old.fy.y*W_t[k][i].fy.y + W_old.fy.z*W_t[k][i].fy.z;
//                tensdottens.z += W_old.fz.x*W_t[k][i].fz.x + W_old.fz.y*W_t[k][i].fz.y + W_old.fz.z*W_t[k][i].fz.z;
//            }
//
//            gamma.x = sigma_old.x/DT2/tensdottens.x;
//            gamma.y = sigma_old.y/DT2/tensdottens.y;
//            gamma.z = sigma_old.z/DT2/tensdottens.z;
//
//            for (i=0; i<_NPART; i++) {
//
//                r_OLD[i].x -= DT2*(gamma.x*W_t[k][i].fx.x + gamma.y*W_t[k][i].fy.x + gamma.z*W_t[k][i].fz.x);
//                r_OLD[i].y -= DT2*(gamma.x*W_t[k][i].fx.y + gamma.y*W_t[k][i].fy.y + gamma.z*W_t[k][i].fz.y);
//                r_OLD[i].z -= DT2*(gamma.x*W_t[k][i].fx.z + gamma.y*W_t[k][i].fy.z + gamma.z*W_t[k][i].fz.z);
//            }
//        }
//
//        count++;
//
//        if (count>_MAX_ITER) {
//
//            printf("SHAKE.c -> ERROR: Iteration time limit exceeded. Convergence not reached!!\n");
//            exit(EXIT_FAILURE);
//
//        }else{
//
//            if (_DEBUG && _D_SHAKE) printf("Iteration: %d\tmaxsigma = (%.4e, %.4e, %.4e)\n", count, maxsigma.x, maxsigma.y, maxsigma.z);
//        }
//
//    } while (maxsigma.x>_TOL || maxsigma.y>_TOL || maxsigma.z>_TOL);
//
//    printf("Convergence of SHAKE reached after %d iteration(s)\n", count);
//    printf("maxsigma = (%.4e, %.4e, %.4e)\n\n", maxsigma.x, maxsigma.y, maxsigma.z);
//}

void TEST(struct point r_t[], struct point r_OLD[]){
    
//    struct tensor TEST = {0};
    
//    TEST.fx.x = 1.;
    
    //    for (k=0; k<_NPART; k++) {
    //
    //        for (i=0; i<_NPART; i++) {
    //
    //            W_t[k][i] = ConstTens(r_t, r_OLD, k, i); // TO BE COMPUTED FOR r(t)
    //        }
    //    }
    
    //    for (i=0; i<_NPART; i++) {
    //
    //        r_OLD[i].x += 2*r_t[i].x;
    //        r_OLD[i].y += 2*r_t[i].y;
    //        r_OLD[i].z += 2*r_t[i].z;
    //    }
    
    printf("TEST\n");
}
