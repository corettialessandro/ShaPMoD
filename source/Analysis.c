//
//  Analysis.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/21/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "Analysis.h"
#include "cells.h"

struct point Velocity(struct point r_tm1, struct point r_tp1){

    struct point v;

    double DT2 = 2.*DT;

    v.x = (r_tp1.x - r_tm1.x)/DT2;
    v.y = (r_tp1.y - r_tm1.y)/DT2;
    v.z = (r_tp1.z - r_tm1.z)/DT2;

    return v;
}

struct point CMVelocity(struct point v[]){

    struct point CMV = {0};

    int i = 0;
    double mi;

    for (i=0; i<NPART; i++){

        mi = M[INDX[i]];

        CMV.x += mi*v[i].x;
        CMV.y += mi*v[i].y;
        CMV.z += mi*v[i].z;
    }

    CMV.x /= MTOT;
    CMV.y /= MTOT;
    CMV.z /= MTOT;

    return CMV;
}

double EnerKin(struct point v[]){

    double EKin = 0;

    int i=0;

    for (i=0; i<NPART; i++) {

        EKin += M[INDX[i]]*modsq(v[i]);
    }

    return EKin*.5;
}

double EnerPot_HM(struct point r[]){

    double EPot = 0;

    int i, j, indx_i, indx_j, indx_int;
    double kQiQj, Aij, Cij, Dij;
    double CC_r, CC_r2, CC_r6, CC_r8;
    struct point CC_d;

    double cut_corr = 0;

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

        for (j=0; j<NPART; j++) {

            if (i!=j) {

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                //Core-Core interactions
//                CC_d = d_rirj(r[i], r[j]);
                CC_d = Distance(r[i], r[j]);
                CC_r = mod(CC_d);

                if (CC_r<LR_RCUT || R_CUT == -1) {

                    kQiQj = _COULOMB*(Q[indx_i] + CHI[indx_i])*(Q[indx_j] + CHI[indx_j]);

                    EPot += kQiQj/CC_r - kQiQj*LR_VCUT - kQiQj*LR_DVCUT*(CC_r - LR_RCUT) - .5*kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)*(CC_r - LR_RCUT);
                }

                if (CC_r<RCUT || R_CUT == -1) {

                    Aij = A[indx_int];
                    Cij = C[indx_int];
                    Dij = D[indx_int];

                    CC_r2 = CC_r*CC_r;
                    CC_r6 = CC_r2*CC_r2*CC_r2;
                    CC_r8 = CC_r6*CC_r2;

                    EPot += Aij*exp(-CC_r/LAMBDA) -Cij/CC_r6 -Dij/CC_r8 - (C_VCUT[indx_int] + S_VCUT[indx_int]) - (C_DVCUT[indx_int] + S_DVCUT[indx_int])*(CC_r - RCUT) - .5*(C_DDVCUT[indx_int] + S_DDVCUT[indx_int])*(CC_r - RCUT)*(CC_r - RCUT);
                }
            }
        }
    }

    for (i=0; i<NINTER; i++) {

        cut_corr += C_VTAIL[i] + S_VTAIL[i];
    }

    return .5*EPot + cut_corr;
}

double EnerPot_Jac(struct point r[], struct point rho[]){

    double EPot = 0;

    int i, j, indx_i, indx_j, indx_int;
    double kQiQj, kChiiChij, kQiChij, kChiiQj, Aij, Cij, Dij;
    double CS_r, SC_r, SC_r2, CC_r, CC_r2, SS_r, CC_r6, CC_r8;
    struct point CS_d, SC_d, CC_d, SS_d;

    double cut_corr = 0;

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

        //Shell-Core self interaction
        SC_d.x = (rho[i].x - r[i].x);
        SC_d.y = (rho[i].y - r[i].y);
        SC_d.z = (rho[i].z - r[i].z);

//        SC_d = d_rhoirj(rho[i], r[i], r[i]);
//        SC_d = Distance(rho[i], r[i]);
        SC_r2 = modsq(SC_d);

        EPot += K[indx_i]*SC_r2;

        if (DEBUG_FLAG && _D_ENERGY) printf("Epot = %lf\n", EPot);

        for (j=0; j<NPART; j++) {

            if (i!=j) {

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                //Core-Core interactions
//                CC_d = d_rirj(r[i], r[j]);
                CC_d = Distance(r[i], r[j]);
                CC_r = mod(CC_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("i = %d\tj = %d\tindx_i = %d\tindx_j = %d\tindx_int = %d\nCC_d = (%lf, %lf, %lf)\tCC_r = %lf\n", i, j, indx_i, indx_j, indx_int, CC_d.x, CC_d.y, CC_d.z, CC_r);

                //Long range C-C interactions
                if (CC_r < LR_RCUT || R_CUT == -1) {

                    kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                    EPot += kQiQj/CC_r - kQiQj*LR_VCUT - kQiQj*LR_DVCUT*(CC_r - LR_RCUT) - .5*kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)*(CC_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("QiQj = %lf\tCC Coulomb term = %lf\n", kQiQj, kQiQj/CC_r);
                }

                //Short range C-C interactions
                if (CC_r < RCUT || R_CUT == -1) {

                    Cij = C[indx_int];
                    Dij = D[indx_int];

                    CC_r2 = CC_r*CC_r;
                    CC_r6 = CC_r2*CC_r2*CC_r2;
                    CC_r8 = CC_r6*CC_r2;

                    EPot += -Cij/CC_r6 -Dij/CC_r8 - C_VCUT[indx_int] - C_DVCUT[indx_int]*(CC_r - RCUT) - .5*C_DDVCUT[indx_int]*(CC_r - RCUT)*(CC_r - RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("Cij = %lf\tDij = %lf\tCC dispersion term = %lf\n", Cij, Dij, -Cij/CC_r6 -Dij/CC_r8);
                }

                //Shell-Shell interactions
//                SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("SS_d = (%lf, %lf, %lf)\tSS_r = %lf\n", SS_d.x, SS_d.y, SS_d.z, SS_r);

                //Long range S-S interactions
                if (SS_r < LR_RCUT || R_CUT == -1){

                    kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                    EPot += kChiiChij/SS_r - kChiiChij*LR_VCUT - kChiiChij*LR_DVCUT*(SS_r - LR_RCUT) - .5*kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)*(SS_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("ChiiChij = %lf\tSS Coulomb term = %lf\n", kChiiChij, kChiiChij/SS_r);
                }

                //Short range S-S interactions
                if (SS_r < RCUT || R_CUT == -1) {

                    Aij = A[indx_int];

                    EPot += Aij*exp(-SS_r/LAMBDA) - S_VCUT[indx_int] - S_DVCUT[indx_int]*(SS_r - RCUT) - .5*S_DDVCUT[indx_int]*(SS_r - RCUT)*(SS_r - RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("Aij = %lf\tSS repulsion term = %lf\n", Aij, Aij*exp(-SS_r/LAMBDA));
                }

                //Core-Shell interactions
//                CS_d = d_rirhoj(r[i], rho[j], r[j]);
                CS_d = Distance(r[i], rho[j]);
                CS_r = mod(CS_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("CS_d = (%lf, %lf, %lf)\tCS_r = %lf\n", CS_d.x, CS_d.y, CS_d.z, CS_r);

                if (CS_r < LR_RCUT || R_CUT == -1){

                    kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                    EPot += kQiChij/CS_r - kQiChij*LR_VCUT - kQiChij*LR_DVCUT*(CS_r - LR_RCUT) - .5*kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)*(CS_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("QiChij = %lf\tCS Coulomb term = %lf\n", kQiChij, kQiChij/CS_r);
                }

                //Shell-Core interactions
//                SC_d = d_rhoirj(rho[i], r[i], r[j]);
                SC_d = Distance(rho[i], r[j]);
                SC_r = mod(SC_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("SC_d = (%lf, %lf, %lf)\tSC_r = %lf\n", SC_d.x, SC_d.y, SC_d.z, SC_r);

                if (SC_r < LR_RCUT || R_CUT == -1){

                    kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                    EPot += kChiiQj/SC_r - kChiiQj*LR_VCUT - kChiiQj*LR_DVCUT*(SC_r - LR_RCUT) - .5*kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)*(SC_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("ChiiQj = %lf\tSC Coulomb term = %lf\n", kChiiQj, kChiiQj/SC_r);
                }
            }

            if (DEBUG_FLAG && _D_ENERGY) printf("EPot = %lf\n", EPot);
        }
    }

    for (i=0; i<NINTER; i++) {

        cut_corr += C_VTAIL[i] + S_VTAIL[i];
    }

    return .5*EPot + cut_corr;
}

double EnerPot_Cicc(struct point r[], struct point rho[]){

    double EPot = 0;

    int i, j, indx_i, indx_j, indx_int;
    //    int n1, n2;
    double kQiQj, kChiiChij, kQiChij, kChiiQj, Aij, Cij, Dij;
    double SC_r, CS_r, SC_r2, CC_r, CC_r2, SS_r, CC_r6, CC_r8;
    struct point SC_d, CS_d, CC_d, SS_d;

    double cut_corr = 0;

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

        SC_d.x = (rho[i].x - r[i].x);
        SC_d.y = (rho[i].y - r[i].y);
        SC_d.z = (rho[i].z - r[i].z);

//        SC_d = d_rhoirj(rho[i], r[i], r[i]);
//        SC_d = Distance(rho[i], r[i]);
        SC_r2 = modsq(SC_d);

        EPot += K[indx_i]*SC_r2;

        if (DEBUG_FLAG && _D_ENERGY) printf("Epot = %lf\n", EPot);

        for (j=0; j<NPART; j++) {

            if (i!=j) {

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                //Core-Core interactions
//                CC_d = d_rirj(r[i], r[j]);
                CC_d = Distance(r[i], r[j]);
                CC_r = mod(CC_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("i = %d\tj = %d\tindx_i = %d\tindx_j = %d\tindx_int = %d\nCC_d = (%lf, %lf, %lf)\tCC_r = %lf\n", i, j, indx_i, indx_j, indx_int, CC_d.x, CC_d.y, CC_d.z, CC_r);

                if (CC_r<LR_RCUT || R_CUT == -1) {

                    kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                    EPot += kQiQj/CC_r - kQiQj*LR_VCUT - kQiQj*LR_DVCUT*(CC_r - LR_RCUT) - .5*kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)*(CC_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("QiQj = %lf\tCC Coulomb term = %lf\n", kQiQj, kQiQj/CC_r);
                }

                if (CC_r<RCUT || R_CUT == -1) {

                    Aij = A[indx_int];
                    Cij = C[indx_int];
                    Dij = D[indx_int];

                    CC_r2 = CC_r*CC_r;
                    CC_r6 = CC_r2*CC_r2*CC_r2;
                    CC_r8 = CC_r6*CC_r2;

                    EPot += Aij*exp(-CC_r/LAMBDA) -Cij/CC_r6 -Dij/CC_r8 - C_VCUT[indx_int] - C_DVCUT[indx_int]*(CC_r - RCUT) - .5*C_DDVCUT[indx_int]*(CC_r - RCUT)*(CC_r - RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("Cij = %lf\tDij = %lf\tCC dispersion term = %lf\n", Cij, Dij, -Cij/CC_r6 -Dij/CC_r8);
                }

                //Shell-Shell interactions
//                SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("SS_d = (%lf, %lf, %lf)\tSS_r = %lf\n", SS_d.x, SS_d.y, SS_d.z, SS_r);

                if (SS_r<LR_RCUT || R_CUT == -1) {

                    kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                    EPot += kChiiChij/SS_r - kChiiChij*LR_VCUT - kChiiChij*LR_DVCUT*(SS_r - LR_RCUT) - .5*kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)*(SS_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("ChiiChij = %lf\tSS Coulomb term = %lf\n", kChiiChij, kChiiChij/SS_r);

                    //                    SS_COUNT[indx_int]++;
                }

                //Core-Shell interactions
//                CS_d = d_rirhoj(r[i], rho[j], r[j]);
                CS_d = Distance(r[i], rho[j]);
                CS_r = mod(CS_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("CS_d = (%lf, %lf, %lf)\tCS_r = %lf\n", CS_d.x, CS_d.y, CS_d.z, CS_r);

                if (CS_r<LR_RCUT || R_CUT == -1) {

                    kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                    EPot += kQiChij/CS_r - kQiChij*LR_VCUT - kQiChij*LR_DVCUT*(CS_r - LR_RCUT) - .5*kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)*(CS_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("QiChij = %lf\tCS Coulomb term = %lf\n", kQiChij, kQiChij/CS_r);

                    //                    CS_COUNT[indx_int]++;
                }

                //Shell-Core interactions
//                SC_d = d_rhoirj(rho[i], r[i], r[j]);
                SC_d = Distance(rho[i], r[j]);
                SC_r = mod(SC_d);

                if (DEBUG_FLAG && _D_ENERGY) printf("SC_d = (%lf, %lf, %lf)\tSC_r = %lf\n", SC_d.x, SC_d.y, SC_d.z, SC_r);

                if (SC_r<LR_RCUT || R_CUT == -1) {

                    kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                    EPot += kChiiQj/SC_r - kChiiQj*LR_VCUT - kChiiQj*LR_DVCUT*(SC_r - LR_RCUT) - .5*kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)*(SC_r - LR_RCUT);

                    if (DEBUG_FLAG && _D_ENERGY) printf("ChiiQj = %lf\tSC Coulomb term = %lf\n", kChiiQj, kChiiQj/SC_r);

                    //                    SC_COUNT[indx_int]++;
                }
            }

            if (DEBUG_FLAG && _D_ENERGY) printf("EPot = %lf\n", EPot);
        }
    }

    for (i=0; i<NINTER; i++) {

        //        n1 = (int) floor(i/2.);
        //        n2 = (int) floor((i+1)/2.);

        cut_corr += C_VTAIL[i];
    }

    return .5*EPot + cut_corr;
}

double EnerPot_WAC(struct point r[]){

    double EPot = 0;

    int i, j, indx_i, indx_j, indx_int;
    //    int n1, n2;
    double SC_r, CS_r, SC_r2, CC_r, CC_r2, SS_r, CC_r3, CC_r6, CC_r8, CC_r12;
    struct point SC_d, CS_d, CC_d, SS_d;
    double LJ_sigma6, LJ_sigma12;

    double cut_corr = 0;

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

//        SC_d = d_rhoirj(rho[i], r[i], r[i]);
//        SC_d = Distance(rho[i], r[i]);

        int p;
        int neighlist[50];
        List_Of_Neighs(i,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {
            j = neighlist[p];
            if (i!=j) {

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
                LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);

                //Core-Core interactions
//              CC_d = d_rirj(r[i], r[j]);
                CC_d = Distance(r[i], r[j]);
                CC_r = mod(CC_d);

                if (CC_r <= LJRCUT[indx_i][indx_j]) {

                    CC_r2 = CC_r*CC_r;
                    CC_r3 = CC_r*CC_r2;
                    CC_r6 = CC_r3*CC_r3;
                    CC_r12 = CC_r6*CC_r6;

                    EPot += 4.*(CC_r12*LJ_sigma12 - CC_r6*LJ_sigma6)*LJEPS[indx_i][indx_j];

                }
            }

        }
    }

    // for (i=0; i<NINTER; i++) {
    //
    //     //        n1 = (int) floor(i/2.);
    //     //        n2 = (int) floor((i+1)/2.);
    //
    //     cut_corr += C_VTAIL[i];
    // }

    return EPot;
}

double EnerTot(struct point r[], struct point rho[], struct point v[]){

    double ETot = 0;

    double EKin, EPot = 0;

    EKin = EnerKin(v);

    if (POT == 'J') {

        EPot = EnerPot_Jac(r, rho);

    }else if (POT == 'C') {

        EPot = EnerPot_Cicc(r, rho);

    }else if (POT == 'C') {

        EPot = EnerPot_WAC(r);

    }

    ETot = EKin + EPot;

    return ETot;
}

double Temperature(struct point v[]){

    double T = 0;

    int i;
    double g = 3*(NPART - 1.);

    for (i=0; i<NPART; i++) {

        T += M[INDX[i]]*modsq(v[i]);
    }

    return T/g;
}

double Pressure_Jac(struct point r[], struct point rho[]){

    double P = 0;

    int i;
    double W_C, W_S;

    for (i=0; i<NPART; i++) {

        W_C = CoreVirial_Jac(r, rho, i);
        W_S = ShellVirial_Jac(rho, r, i);

        P += W_C + W_S;
    }

    return P/6./LBOX/LBOX/LBOX;
}

double Pressure_Cicc(struct point r[], struct point rho[]){

    double P = 0;

    int i;
    double W_C, W_S;

    for (i=0; i<NPART; i++) {

        W_C = CoreVirial_Cicc(r, rho, i);
        W_S = ShellVirial_Cicc(rho, r, i);

        P += W_C + W_S;
    }

    return P/6./LBOX/LBOX/LBOX;
}

void GofR(struct point r[]){

    int i, j, indx_i, indx_j, indx_int, indx_g;

    struct point CC_d;
    double CC_r;

    for (i=0; i<NPART-1; i++) {

        indx_i = INDX[i];

        for (j=i+1; j<NPART; j++) {

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

//            CC_d = d_rirj(r[i], r[j]);
            CC_d = Distance(r[i], r[j]);
            CC_r = mod(CC_d);

            if (CC_r<.5*LBOX) {

                indx_g = (int) floor(CC_r/DG);
                GCOUNT[indx_int][indx_g] += 2;
            }
        }
    }

    NGCOUNT++;
}

struct point Angular_Momentum(struct point r[], struct point v[]){

    struct point Omega_tot = {0}, Omega;

    struct point pi;

    int i;

    for(i=0; i<NPART; i++){

        pi.x = M[INDX[i]]*v[i].x;
        pi.y = M[INDX[i]]*v[i].y;
        pi.z = M[INDX[i]]*v[i].z;

        Omega = Vector(r[i], pi);

        Omega_tot.x += Omega.x;
        Omega_tot.y += Omega.y;
        Omega_tot.z += Omega.z;
    }

    return Omega_tot;
}

void Analyse(int timestep, struct point r[], struct point rho[], struct point v[], char therm){

    int i;
    double ETot, EKin, EPot = 0, EPol, Temp, Press = 0;
    struct point CMV, AnMom;

    EKin = EnerKin(v)*_E_CONV;
    ITEMP = Temp = Temperature(v);

    for (i=0; i<NINTER; i++) {

        Press = C_PTAIL[i] + S_PTAIL[i];
    }

    if (POT == 'J') {

        EPot = EnerPot_Jac(r, rho)*_E_CONV;
        Press += (DENSITY*_KB*Temp + Pressure_Jac(r, rho))*_PRESS_CONV;

        //printf("P_IG = %lf\tP_VIR = %lf\tP_TOT = %lf\n", DENSITY*_KB*Temp, Pressure_Jac(r, rho), DENSITY*_KB*Temp + Pressure_Jac(r, rho));

    } else if (POT == 'C'){

        EPot = EnerPot_Cicc(r, rho)*_E_CONV;
        Press += (DENSITY*_KB*Temp + Pressure_Cicc(r, rho))*_PRESS_CONV;

    } else if (POT == 'W'){

        EPot = EnerPot_WAC(r)*_E_CONV;
        Press += 0;
    }

    EPol = EPot - EnerPot_HM(r)*_E_CONV;

    Temp *= _TEMP_CONV;
    ETot = EKin + EPot;
    CMV = CMVelocity(v);
    AnMom = Angular_Momentum(r, v);

    NANCOUNT++;

    EKINAVG = OTRAvg(NANCOUNT, EKin, EKINAVG);
    EPOTAVG = OTRAvg(NANCOUNT, EPot, EPOTAVG);
    ETOTAVG = EKINAVG + EPOTAVG;
    EPOLAVG = OTRAvg(NANCOUNT, EPol, EPOLAVG);
    TEMPAVG = OTRAvg(NANCOUNT, Temp, TEMPAVG);
    PRESSAVG = OTRAvg(NANCOUNT, Press, PRESSAVG);
    ANMOMAVG.x = OTRAvg(NANCOUNT, AnMom.x, ANMOMAVG.x);
    ANMOMAVG.y = OTRAvg(NANCOUNT, AnMom.y, ANMOMAVG.y);
    ANMOMAVG.z = OTRAvg(NANCOUNT, AnMom.z, ANMOMAVG.z);

    ETOTAVGSQ = OTRAvg(NANCOUNT, ETot*ETot, ETOTAVGSQ);
    EKINAVGSQ = OTRAvg(NANCOUNT, EKin*EKin, EKINAVGSQ);
    EPOTAVGSQ = OTRAvg(NANCOUNT, EPot*EPot, EPOTAVGSQ);
    EPOLAVGSQ = OTRAvg(NANCOUNT, EPol*EPol, EPOLAVGSQ);
    TEMPAVGSQ = OTRAvg(NANCOUNT, Temp*Temp, TEMPAVGSQ);
    PRESSAVGSQ = OTRAvg(NANCOUNT, Press*Press, PRESSAVGSQ);

    if (timestep % IANSCREEN == 0) printf("%d\t%.12e\t%.7e\t%.7e\t%.7e\t%c%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.7e\t%.7e\t%.7e\t%d\t%.4e\n", timestep, ETot, EKin, EPot, 100.*EPol/EPot, therm, Temp, Press, CMV.x, CMV.y, CMV.z,  AnMom.x, AnMom.y, AnMom.z, SR_ITERS, SR_DISCR);

    Write_Analysis(timestep, ETot, EKin, EPot, EPol, therm, Temp, Press, CMV, AnMom, SR_ITERS, SR_DISCR);
}
