//
//  Forces.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "Forces.h"
#include "cells.h"

struct point CoreForce_Jac(struct point r[], struct point rho[], int i){

    struct point F;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kQiQj, kQiChij, Cij, Dij;
    struct point CC_d, CS_d;
    double CC_r, CC_r2, CC_r3, CC_r8, CC_r10;
    double CS_r, CS_r3;

    CS_d.x = (r[i].x - rho[i].x);
    CS_d.y = (r[i].y - rho[i].y);
    CS_d.z = (r[i].z - rho[i].z);

//    CS_d = d_rirhoj(r[i], rho[i], r[i]);
//    CS_d = Distance(r[i], rho[i]);

    if (DEBUG_FLAG && _D_DIST) printf("r[%d] = (%.4e, %.4e, %.4e), rho[%d] = (%.4e, %.4e, %.4e), CS_d = (%.4e, %.4e, %.4e)\n", i, r[i].x, r[i].y, r[i].z, i, rho[i].x, rho[i].y, rho[i].z, CS_d.x, CS_d.y, CS_d.z);

    F.x = -ki*CS_d.x;
    F.y = -ki*CS_d.y;
    F.z = -ki*CS_d.z;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Core-Core interactions
//            CC_d = d_rirj(r[i], r[j]);
            CC_d = Distance(r[i], r[j]);
            CC_r = mod(CC_d);

            if (DEBUG_FLAG && _D_DIST) printf("CC_d[%d][%d] = (%.4e, %.4e, %.4e)\nCC_r[%d][%d] = %.4e\n", i, j, CC_d.x, CC_d.y, CC_d.z, i, j, CC_r);

            if (CC_r < LR_RCUT || R_CUT == -1){

                CC_r2 = CC_r*CC_r;
                CC_r3 = CC_r*CC_r2;

                kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                F.x += kQiQj/CC_r3*CC_d.x + kQiQj*LR_DVCUT/CC_r*CC_d.x + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.x; //coulomb part
                F.y += kQiQj/CC_r3*CC_d.y + kQiQj*LR_DVCUT/CC_r*CC_d.y + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.y;
                F.z += kQiQj/CC_r3*CC_d.z + kQiQj*LR_DVCUT/CC_r*CC_d.z + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.z;
            }

            if (CC_r < RCUT || R_CUT == -1) {

                CC_r2 = CC_r*CC_r;

                Cij = C[indx_int];
                Dij = D[indx_int];

                CC_r8 = CC_r2*CC_r2*CC_r2*CC_r2;
                CC_r10 = CC_r8*CC_r2;

                F.x += - 6.*Cij/CC_r8*CC_d.x - 8.*Dij/CC_r10*CC_d.x + C_DVCUT[indx_int]/CC_r*CC_d.x + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.x; //dispersion terms
                F.y += - 6.*Cij/CC_r8*CC_d.y - 8.*Dij/CC_r10*CC_d.y + C_DVCUT[indx_int]/CC_r*CC_d.y + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.y;
                F.z += - 6.*Cij/CC_r8*CC_d.z - 8.*Dij/CC_r10*CC_d.z + C_DVCUT[indx_int]/CC_r*CC_d.z + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.z;
            }

            //Core-Shell interactions
//            CS_d = d_rirhoj(r[i], rho[j], r[j]);
            CS_d = Distance(r[i], rho[j]);
            CS_r = mod(CS_d);

            if (DEBUG_FLAG && _D_DIST) printf("CS_d[%d][%d] = (%.4e, %.4e, %.4e)\nCS_r[%d][%d] = %.4e\n", i, j, CS_d.x, CS_d.y, CS_d.z, i, j, CS_r);

            if (CS_r < LR_RCUT || R_CUT == -1){

                kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                CS_r3 = CS_r*CS_r*CS_r;

                F.x += kQiChij/CS_r3*CS_d.x + kQiChij*LR_DVCUT/CS_r*CS_d.x + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.x;
                F.y += kQiChij/CS_r3*CS_d.y + kQiChij*LR_DVCUT/CS_r*CS_d.y + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.y;
                F.z += kQiChij/CS_r3*CS_d.z + kQiChij*LR_DVCUT/CS_r*CS_d.z + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.z;
            }
        }
    }

    return F;
}

struct point ShellForce_Jac(struct point rho[], struct point r[], int i){

    struct point F;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj, Aij;
    struct point SS_d, SC_d;
    double SS_r, SS_r3;
    double SC_r, SC_r3;

    SC_d.x = (rho[i].x - r[i].x);
    SC_d.y = (rho[i].y - r[i].y);
    SC_d.z = (rho[i].z - r[i].z);

//    SC_d = d_rhoirj(rho[i], r[i], r[i]);
//    SC_d = Distance(rho[i], r[i]);

    if (DEBUG_FLAG && _D_DIST) printf("rho[%d] = (%.4e, %.4e, %.4e), r[%d] = (%.4e, %.4e, %.4e), SC_d = (%.4e, %.4e, %.4e)\n", i, rho[i].x, rho[i].y, rho[i].z, i, r[i].x, r[i].y, r[i].z, SC_d.x, SC_d.y, SC_d.z);

    F.x = -ki*SC_d.x; //forces harmonic term
    F.y = -ki*SC_d.y;
    F.z = -ki*SC_d.z;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Shell-Shell interactions
//            SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
            SS_d = Distance(rho[i], rho[j]);
            SS_r = mod(SS_d);

            if (DEBUG_FLAG && _D_DIST) printf("SS_d[%d][%d] = (%.4e, %.4e, %.4e)\nSS_r[%d][%d] = %.4e\n", i, j, SS_d.x, SS_d.y, SS_d.z, i, j, SS_r);

            if (SS_r < LR_RCUT || R_CUT == -1){

                kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                SS_r3 = SS_r*SS_r*SS_r;

                F.x += kChiiChij/SS_r3*SS_d.x + kChiiChij*LR_DVCUT/SS_r*SS_d.x + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.x; //coulomb part
                F.y += kChiiChij/SS_r3*SS_d.y + kChiiChij*LR_DVCUT/SS_r*SS_d.y + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.y;
                F.z += kChiiChij/SS_r3*SS_d.z + kChiiChij*LR_DVCUT/SS_r*SS_d.z + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.z;
            }

            if (SS_r < RCUT || R_CUT == -1) {

                Aij = A[indx_int];

                F.x += Aij/LAMBDA*SS_d.x/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.x + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.x; //exponential part
                F.y += Aij/LAMBDA*SS_d.y/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.y + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.y;
                F.z += Aij/LAMBDA*SS_d.z/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.z + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.z;

                //                printf("F_rep = (%lf, %lf, %lf)\n", Aij/LAMBDA*SS_d.x/SS_r*exp(-SS_r/LAMBDA), Aij/LAMBDA*SS_d.y/SS_r*exp(-SS_r/LAMBDA), Aij/LAMBDA*SS_d.z/SS_r*exp(-SS_r/LAMBDA));
            }

            //Shell-Core interactions
//            SC_d = d_rhoirj(rho[i], r[i], r[j]);
            SC_d = Distance(rho[i], r[j]);
            SC_r = mod(SC_d);

            if (DEBUG_FLAG && _D_DIST) printf("SC_d[%d][%d] = (%.4e, %.4e, %.4e)\nSC_r[%d][%d] = %.4e\n", i, j, SC_d.x, SC_d.y, SC_d.z, i, j, SC_r);

            if (SC_r < LR_RCUT || R_CUT == -1){

                kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                SC_r3 = SC_r*SC_r*SC_r;

                F.x += kChiiQj/SC_r3*SC_d.x + kChiiQj*LR_DVCUT/SC_r*SC_d.x + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.x;
                F.y += kChiiQj/SC_r3*SC_d.y + kChiiQj*LR_DVCUT/SC_r*SC_d.y + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.y;
                F.z += kChiiQj/SC_r3*SC_d.z + kChiiQj*LR_DVCUT/SC_r*SC_d.z + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.z;
            }
        }

        if (DEBUG_FLAG && _D_FORCES) printf("+F[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, F.x, F.y, F.z);
    }

    return F;
}

struct point CoreForce_Cicc(struct point r[], struct point rho[], int i){

    struct point F;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kQiQj, kQiChij, Aij, Cij, Dij;
    struct point CC_d, CS_d;
    double CC_r, CC_r2, CC_r3, CC_r8, CC_r10;
    double CS_r, CS_r3;

    CS_d.x = (r[i].x - rho[i].x);
    CS_d.y = (r[i].y - rho[i].y);
    CS_d.z = (r[i].z - rho[i].z);

//    CS_d = Distance(r[i], rho[i]);
//    CS_d = d_rirhoj(r[i], rho[i], r[i]);

    if (DEBUG_FLAG && _D_DIST) {
        printf("r[%d] = (%.4e, %.4e, %.4e), rho[%d] = (%.4e, %.4e, %.4e), CS_d = (%.4e, %.4e, %.4e)\n", i, r[i].x, r[i].y, r[i].z, i, rho[i].x, rho[i].y, rho[i].z, CS_d.x, CS_d.y, CS_d.z);
    }

    F.x = -ki*CS_d.x;
    F.y = -ki*CS_d.y;
    F.z = -ki*CS_d.z;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Core-Core interactions
//            CC_d = d_rirj(r[i], r[j]);
            CC_d = Distance(r[i], r[j]);
            if (DEBUG_FLAG && _D_DIST) printf("CC_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, CC_d.x, CC_d.y, CC_d.z);

            CC_r = mod(CC_d);
            if (DEBUG_FLAG && _D_DIST) printf("CC_r[%d][%d] = %.4e\n", i, j, CC_r);

            if (CC_r<LR_RCUT || R_CUT == -1) {

                kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                CC_r2 = CC_r*CC_r;
                CC_r3 = CC_r*CC_r2;

                F.x += kQiQj/CC_r3*CC_d.x + kQiQj*LR_DVCUT/CC_r*CC_d.x + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.x;
                F.y += kQiQj/CC_r3*CC_d.y + kQiQj*LR_DVCUT/CC_r*CC_d.y + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.y;
                F.z += kQiQj/CC_r3*CC_d.z + kQiQj*LR_DVCUT/CC_r*CC_d.z + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.z;
            }

            if (CC_r<RCUT || R_CUT == -1) {

                CC_r2 = CC_r*CC_r;

                Aij = A[indx_int];
                Cij = C[indx_int];
                Dij = D[indx_int];

                CC_r8 = CC_r2*CC_r2*CC_r2*CC_r2;
                CC_r10 = CC_r8*CC_r2;

                F.x += Aij/LAMBDA*CC_d.x/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.x - 8.*Dij/CC_r10*CC_d.x + C_DVCUT[indx_int]/CC_r*CC_d.x + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.x;
                F.y += Aij/LAMBDA*CC_d.y/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.y - 8.*Dij/CC_r10*CC_d.y + C_DVCUT[indx_int]/CC_r*CC_d.y + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.y;
                F.z += Aij/LAMBDA*CC_d.z/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.z - 8.*Dij/CC_r10*CC_d.z + C_DVCUT[indx_int]/CC_r*CC_d.z + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.z;
            }

            //Core-Shell interactions
//            CS_d = d_rirhoj(r[i], rho[j], r[j]);
            CS_d = Distance(r[i], rho[j]);
            if (DEBUG_FLAG && _D_DIST) printf("CS_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, CS_d.x, CS_d.y, CS_d.z);

            CS_r = mod(CS_d);
            if (DEBUG_FLAG && _D_DIST) printf("CS_r[%d][%d] = %.4e\n", i, j, CS_r);

            if (CS_r<LR_RCUT || R_CUT == -1) {

                kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                CS_r3 = CS_r*CS_r*CS_r;

                F.x += kQiChij/CS_r3*CS_d.x + kQiChij*LR_DVCUT/CS_r*CS_d.x + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.x;
                F.y += kQiChij/CS_r3*CS_d.y + kQiChij*LR_DVCUT/CS_r*CS_d.y + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.y;
                F.z += kQiChij/CS_r3*CS_d.z + kQiChij*LR_DVCUT/CS_r*CS_d.z + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.z;
            }
        }
    }

    return F;
}

struct point ShellForce_Cicc(struct point rho[], struct point r[], int i){

    struct point F;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj;
    struct point SS_d, SC_d;
    double SS_r, SS_r3;
    double SC_r, SC_r3;

    SC_d.x = (rho[i].x - r[i].x);
    SC_d.y = (rho[i].y - r[i].y);
    SC_d.z = (rho[i].z - r[i].z);

//    SC_d = Distance(rho[i], r[i]);
//    SC_d = d_rhoirj(rho[i], r[i], r[i]);

    if (DEBUG_FLAG && _D_DIST){

        printf("rho[%d] = (%.4e, %.4e, %.4e), r[%d] = (%.4e, %.4e, %.4e), SC_d = (%.4e, %.4e, %.4e)\n", i, rho[i].x, rho[i].y, rho[i].z, i, r[i].x, r[i].y, r[i].z, SC_d.x, SC_d.y, SC_d.z);
    }

    F.x = -ki*SC_d.x;
    F.y = -ki*SC_d.y;
    F.z = -ki*SC_d.z;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Shell-Shell interactions
//            SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
            SS_d = Distance(rho[i], rho[j]);
            if (DEBUG_FLAG && _D_DIST) printf("SS_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, SS_d.x, SS_d.y, SS_d.z);

            SS_r = mod(SS_d);
            if (DEBUG_FLAG && _D_DIST) printf("SS_r[%d][%d] = %.4e\n", i, j, SS_r);

            if (SS_r<LR_RCUT || R_CUT == -1) {

                kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                SS_r3 = SS_r*SS_r*SS_r;

                F.x += kChiiChij/SS_r3*SS_d.x + kChiiChij*LR_DVCUT/SS_r*SS_d.x + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.x;
                F.y += kChiiChij/SS_r3*SS_d.y + kChiiChij*LR_DVCUT/SS_r*SS_d.y + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.y;
                F.z += kChiiChij/SS_r3*SS_d.z + kChiiChij*LR_DVCUT/SS_r*SS_d.z + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.z;
            }

            //Shell-Core interactions
//            SC_d = d_rhoirj(rho[i], r[i], r[j]);
            SC_d = Distance(rho[i], r[j]);
            if (DEBUG_FLAG && _D_DIST) printf("SC_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, SC_d.x, SC_d.y, SC_d.z);

            SC_r = mod(SC_d);
            if (DEBUG_FLAG && _D_DIST) printf("SC_r[%d][%d] = %.4e\n", i, j, SC_r);

            if (SC_r<LR_RCUT || R_CUT == -1) {

                kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                SC_r3 = SC_r*SC_r*SC_r;

                F.x += kChiiQj/SC_r3*SC_d.x + kChiiQj*LR_DVCUT/SC_r*SC_d.x + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.x;
                F.y += kChiiQj/SC_r3*SC_d.y + kChiiQj*LR_DVCUT/SC_r*SC_d.y + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.y;
                F.z += kChiiQj/SC_r3*SC_d.z + kChiiQj*LR_DVCUT/SC_r*SC_d.z + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.z;
            }
        }

        if (DEBUG_FLAG && _D_FORCES) printf("+F[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, F.x, F.y, F.z);
    }

    return F;
}

struct point Force_WCA(struct point r[], int i) {

    struct point F;

    double rCUT;
    double lround;
    double Rround, R2round;
    int j, indx_i = INDX[i], indx_j, indx_int;
    struct point CC_d;
    double CC_r, r2inv, r6inv, rc2inv, rc6inv, ljatrc, forcelj, fpair, ljrc;
    double LJ_sigma6, LJ_sigma12;
    int nn;

    //CS_d = d_rirhoj(r[i], rho[i], r[i]);
    //CS_d = Distance(r[i], rho[i]);
    F.x = 0.;
    F.y = 0.;
    F.z = 0.;

    int p;
    int neighlist[1000];
    List_Of_Neighs(i,neighlist,1);
    nn = 0;
    //if (i == 752) printf("Const %d\n", neighlist[0]);
    //if (i==100) printf("Neighs (cells):");
    //for (j=0; j<NPART; j++){
    for (p=1;p<=neighlist[0];p++) {
        j = neighlist[p];
        if (i!=j) {
          indx_j = INDX[j];
          indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

          rCUT = LJRCUT[indx_i][indx_j];
          lround = LJLROUND[indx_i][indx_j];
          LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
          LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);

          //Core-Core interactions
          //CC_d = d_rirj(r[i], r[j]);
          CC_d = Distance(r[i], r[j]);
          CC_r = mod(CC_d);
          //if (i == 752) printf("Force : r = %.4e and rcut = %.4e\n", CC_r, rCUT);
          // printf("CC_r = %.4e \n", CC_r);
          // printf("rCUT = %.4e \n", rCUT);
          // exit(0);

          if (CC_r <= rCUT){
              nn++;
              //if (i==100) printf(" %d",j);

              r2inv = 1.0/(CC_r*CC_r);
              r6inv = r2inv*r2inv*r2inv;
              rc2inv = 1.0/(rCUT*rCUT);
              rc6inv = rc2inv*rc2inv*rc2inv;
              ljatrc = rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
              //ljatrc = 4. * LJEPS[indx_i][indx_j] * rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
              forcelj = 24.0 * LJEPS[indx_i][indx_j]/CC_r * r6inv * (2.0*LJ_sigma12*r6inv - LJ_sigma6);

              //rounding the force

              if (CC_r > (rCUT-lround)) {
                  Rround = (CC_r - rCUT + lround)/lround;
                  //forcelj = 24.0
                  //          * (Rround * (Rround-1.0)/lround * (LJEPS[indx_i][indx_j] * r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + 0.25*ljatrc) //0.25 missing?
                  //          + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * LJEPS[indx_i][indx_j] * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
                  forcelj = 24.0 * LJEPS[indx_i][indx_j]
                            * (Rround * (Rround-1.0)/lround * (r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + ljatrc)
                               + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
              }

              fpair = forcelj/CC_r;

              F.x += CC_d.x*fpair;
              F.y += CC_d.y*fpair;
              F.z += CC_d.z*fpair;

          }

        }
    }
    //if (i==100) printf(" (%d)\n",nn);
/*
    nn=0;
    double fx=F.x, fy=F.y, fz=F.z;
    //if (i==100) printf("Neighs:");
    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            rCUT = LJRCUT[indx_i][indx_j];
            LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
            LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);

            //Core-Core interactions
            //CC_d = d_rirj(r[i], r[j]);
            CC_d = Distance(r[i], r[j]);
            CC_r = mod(CC_d);
            // printf("CC_r = %.4e \n", CC_r);
            // printf("rCUT = %.4e \n", rCUT);
            // exit(0);

            if (CC_r <= rCUT){
                nn++;
                //if (i==100) printf(" %d",j);

                r2inv = 1.0/(CC_r*CC_r);
                r6inv = r2inv*r2inv*r2inv;
                rc2inv = 1.0/(rCUT*rCUT);
                rc6inv = rc2inv*rc2inv*rc2inv;
                ljatrc = rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
                forcelj = 24.0 * LJEPS[indx_i][indx_j]/CC_r * r6inv * (2.0*LJ_sigma12*r6inv - LJ_sigma6);

                //rounding the force

                if (CC_r > (rCUT-lround)) {
                    Rround = (CC_r - rCUT + lround)/lround;
                    forcelj = 24.0 * LJEPS[indx_i][indx_j]
                            * (Rround * (Rround-1.0)/lround * (r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + ljatrc)
                               + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
                }

                fpair = forcelj/CC_r;

                F.x += CC_d.x*fpair;
                F.y += CC_d.y*fpair;
                F.z += CC_d.z*fpair;
                //fx -= CC_d.x*fpair;
                //fy -= CC_d.y*fpair;
                //fz -= CC_d.z*fpair;

            }

        }
    }
*/
    //if (fabs(fx>0.000001)) printf("Interaction %d (%lf,%lf,%lf) along x wrong!\n",i,r[i].x+0.5*LBOX,r[i].y+0.5*LBOX,r[i].z+0.5*LBOX);
    //if (fabs(fy>0.000001)) printf("Interaction %d (%lf,%lf,%lf) along y wrong!\n",i,r[i].x+0.5*LBOX,r[i].y+0.5*LBOX,r[i].z+0.5*LBOX);
    //if (fabs(fz>0.000001)) printf("Interaction %d (%lf,%lf,%lf) along z wrong!\n",i,r[i].x+0.5*LBOX,r[i].y+0.5*LBOX,r[i].z+0.5*LBOX);
    //if (i==100) printf(" (%d)\n",nn);
    //if (i==100) {printf("Force: %lf %lf %lf (%lf %lf %lf)\n",F.x,F.y,F.z,fx,fy,fz);}

    return F;
}

struct point Force_LJ(struct point r[], int i) {

    struct point F;

    double rCUT;
    double lround;
    double Rround, R2round;
    int j, indx_i = INDX[i], indx_j, indx_int;
    struct point CC_d;
    double CC_r, r2inv, r6inv, rc2inv, rc6inv, ljatrc, forcelj, fpair, ljrc;
    double LJ_sigma6, LJ_sigma12;
    int nn;

    //CS_d = d_rirhoj(r[i], rho[i], r[i]);
    //CS_d = Distance(r[i], rho[i]);
    F.x = 0.;
    F.y = 0.;
    F.z = 0.;

    // int p;
    // int neighlist[1000];
    // List_Of_Neighs(i,neighlist,1);
    // nn = 0;
    //if (i == 752) printf("Const %d\n", neighlist[0]);
    //if (i==100) printf("Neighs (cells):");
    for (j=0; j<NPART; j++){
    // for (p=1;p<=neighlist[0];p++) {
    //     j = neighlist[p];
        if (i!=j) {
          indx_j = INDX[j];
          indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

          rCUT = LJRCUT[indx_i][indx_j];
          lround = LJLROUND[indx_i][indx_j];
          LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
          LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);

          //Core-Core interactions
          //CC_d = d_rirj(r[i], r[j]);
          CC_d = Distance(r[i], r[j]);
          CC_r = mod(CC_d);
          //if (i == 752) printf("Force : r = %.4e and rcut = %.4e\n", CC_r, rCUT);
          // printf("CC_r = %.4e \n", CC_r);
          // printf("rCUT = %.4e \n", rCUT);
          // exit(0);


              //if (i==100) printf(" %d",j);

          r2inv = 1.0/(CC_r*CC_r);
          r6inv = r2inv*r2inv*r2inv;
          rc2inv = 1.0/(rCUT*rCUT);
          rc6inv = rc2inv*rc2inv*rc2inv;
          ljatrc = rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
          //ljatrc = 4. * LJEPS[indx_i][indx_j] * rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
          forcelj = 24.0 * LJEPS[indx_i][indx_j]/CC_r * r6inv * (2.0*LJ_sigma12*r6inv - LJ_sigma6);


          fpair = forcelj/CC_r;

          F.x += CC_d.x*fpair;
          F.y += CC_d.y*fpair;
          F.z += CC_d.z*fpair;

        }
    }

    return F;
}

double CoreVirial_Jac(struct point r[], struct point rho[], int i){

    struct point W;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kQiQj, kQiChij, Cij, Dij;
    struct point CC_d, CS_d;
    double CC_r, CC_r2, CC_r3, CC_r8, CC_r10;
    double CS_r, CS_r3;

    CS_d.x = (r[i].x - rho[i].x);
    CS_d.y = (r[i].y - rho[i].y);
    CS_d.z = (r[i].z - rho[i].z);

    //    CS_d = d_rirhoj(r[i], rho[i], r[i]);
    //    CS_d = Distance(r[i], rho[i]);

    if (DEBUG_FLAG && _D_DIST) printf("r[%d] = (%.4e, %.4e, %.4e), rho[%d] = (%.4e, %.4e, %.4e), CS_d = (%.4e, %.4e, %.4e)\n", i, r[i].x, r[i].y, r[i].z, i, rho[i].x, rho[i].y, rho[i].z, CS_d.x, CS_d.y, CS_d.z);

        W.x = (-ki*CS_d.x)*CS_d.x;
        W.y = (-ki*CS_d.y)*CS_d.y;
        W.z = (-ki*CS_d.z)*CS_d.z;

        for(j=0; j<NPART; j++){

            if (i!=j){

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                //Core-Core interactions
                //            CC_d = d_rirj(r[i], r[j]);
                CC_d = Distance(r[i], r[j]);
                CC_r = mod(CC_d);

                if (DEBUG_FLAG && _D_DIST) printf("CC_d[%d][%d] = (%.4e, %.4e, %.4e)\nCC_r[%d][%d] = %.4e\n", i, j, CC_d.x, CC_d.y, CC_d.z, i, j, CC_r);

                if (CC_r < LR_RCUT || R_CUT == -1){

                    CC_r2 = CC_r*CC_r;
                    CC_r3 = CC_r*CC_r2;

                    kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                    W.x += (kQiQj/CC_r3*CC_d.x + kQiQj*LR_DVCUT/CC_r*CC_d.x + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.x)*CC_d.x;
                    W.y += (kQiQj/CC_r3*CC_d.y + kQiQj*LR_DVCUT/CC_r*CC_d.y + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.y)*CC_d.y;
                    W.z += (kQiQj/CC_r3*CC_d.z + kQiQj*LR_DVCUT/CC_r*CC_d.z + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.z)*CC_d.z;
                }

                if (CC_r < RCUT || R_CUT == -1) {

                    CC_r2 = CC_r*CC_r;

                    Cij = C[indx_int];
                    Dij = D[indx_int];

                    CC_r8 = CC_r2*CC_r2*CC_r2*CC_r2;
                    CC_r10 = CC_r8*CC_r2;

                    W.x += (- 6.*Cij/CC_r8*CC_d.x - 8.*Dij/CC_r10*CC_d.x + C_DVCUT[indx_int]/CC_r*CC_d.x + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.x)*CC_d.x;
                    W.y += (- 6.*Cij/CC_r8*CC_d.y - 8.*Dij/CC_r10*CC_d.y + C_DVCUT[indx_int]/CC_r*CC_d.y + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.y)*CC_d.y;
                    W.z += (- 6.*Cij/CC_r8*CC_d.z - 8.*Dij/CC_r10*CC_d.z + C_DVCUT[indx_int]/CC_r*CC_d.z + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.z)*CC_d.z;
                }

                //Core-Shell interactions
                //            CS_d = d_rirhoj(r[i], rho[j], r[j]);
                CS_d = Distance(r[i], rho[j]);
                CS_r = mod(CS_d);

                if (DEBUG_FLAG && _D_DIST) printf("CS_d[%d][%d] = (%.4e, %.4e, %.4e)\nCS_r[%d][%d] = %.4e\n", i, j, CS_d.x, CS_d.y, CS_d.z, i, j, CS_r);

                if (CS_r < LR_RCUT || R_CUT == -1){

                    kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                    CS_r3 = CS_r*CS_r*CS_r;

                    W.x += (kQiChij/CS_r3*CS_d.x + kQiChij*LR_DVCUT/CS_r*CS_d.x + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.x)*CS_d.x;
                    W.y += (kQiChij/CS_r3*CS_d.y + kQiChij*LR_DVCUT/CS_r*CS_d.y + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.y)*CS_d.y;
                    W.z += (kQiChij/CS_r3*CS_d.z + kQiChij*LR_DVCUT/CS_r*CS_d.z + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.z)*CS_d.z;
                }
            }
        }

    return W.x + W.y + W.z;
}

double ShellVirial_Jac(struct point rho[], struct point r[], int i){

    struct point W;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj, Aij;
    struct point SS_d, SC_d;
    double SS_r, SS_r3;
    double SC_r, SC_r3;

    SC_d.x = (rho[i].x - r[i].x);
    SC_d.y = (rho[i].y - r[i].y);
    SC_d.z = (rho[i].z - r[i].z);

    //    SC_d = d_rhoirj(rho[i], r[i], r[i]);
    //    SC_d = Distance(rho[i], r[i]);

    if (DEBUG_FLAG && _D_DIST) printf("rho[%d] = (%.4e, %.4e, %.4e), r[%d] = (%.4e, %.4e, %.4e), SC_d = (%.4e, %.4e, %.4e)\n", i, rho[i].x, rho[i].y, rho[i].z, i, r[i].x, r[i].y, r[i].z, SC_d.x, SC_d.y, SC_d.z);

        W.x = (-ki*SC_d.x)*SC_d.x;
        W.y = (-ki*SC_d.y)*SC_d.y;
        W.z = (-ki*SC_d.z)*SC_d.z;

        for(j=0; j<NPART; j++){

            if (i!=j){

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                //Shell-Shell interactions
                //            SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (DEBUG_FLAG && _D_DIST) printf("SS_d[%d][%d] = (%.4e, %.4e, %.4e)\nSS_r[%d][%d] = %.4e\n", i, j, SS_d.x, SS_d.y, SS_d.z, i, j, SS_r);

                if (SS_r < LR_RCUT || R_CUT == -1){

                    kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                    SS_r3 = SS_r*SS_r*SS_r;

                    W.x += (kChiiChij/SS_r3*SS_d.x + kChiiChij*LR_DVCUT/SS_r*SS_d.x + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.x)*SS_d.x;
                    W.y += (kChiiChij/SS_r3*SS_d.y + kChiiChij*LR_DVCUT/SS_r*SS_d.y + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.y)*SS_d.y;
                    W.z += (kChiiChij/SS_r3*SS_d.z + kChiiChij*LR_DVCUT/SS_r*SS_d.z + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.z)*SS_d.z;
                }

                if (SS_r < RCUT || R_CUT == -1) {

                    Aij = A[indx_int];

                    W.x += (Aij/LAMBDA*SS_d.x/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.x + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.x)*SS_d.x;
                    W.y += (Aij/LAMBDA*SS_d.y/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.y + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.y)*SS_d.y;
                    W.z += (Aij/LAMBDA*SS_d.z/SS_r*exp(-SS_r/LAMBDA) + S_DVCUT[indx_int]/SS_r*SS_d.z + S_DDVCUT[indx_int]*(SS_r - RCUT)/SS_r*SS_d.z)*SS_d.z;

                    //                printf("F_rep = (%lf, %lf, %lf)\n", Aij/LAMBDA*SS_d.x/SS_r*exp(-SS_r/LAMBDA), Aij/LAMBDA*SS_d.y/SS_r*exp(-SS_r/LAMBDA), Aij/LAMBDA*SS_d.z/SS_r*exp(-SS_r/LAMBDA));
                }

                //Shell-Core interactions
                //            SC_d = d_rhoirj(rho[i], r[i], r[j]);
                SC_d = Distance(rho[i], r[j]);
                SC_r = mod(SC_d);

                if (DEBUG_FLAG && _D_DIST) printf("SC_d[%d][%d] = (%.4e, %.4e, %.4e)\nSC_r[%d][%d] = %.4e\n", i, j, SC_d.x, SC_d.y, SC_d.z, i, j, SC_r);

                if (SC_r < LR_RCUT || R_CUT == -1){

                    kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                    SC_r3 = SC_r*SC_r*SC_r;

                    W.x += (kChiiQj/SC_r3*SC_d.x + kChiiQj*LR_DVCUT/SC_r*SC_d.x + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.x)*SC_d.x;
                    W.y += (kChiiQj/SC_r3*SC_d.y + kChiiQj*LR_DVCUT/SC_r*SC_d.y + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.y)*SS_d.y;
                    W.z += (kChiiQj/SC_r3*SC_d.z + kChiiQj*LR_DVCUT/SC_r*SC_d.z + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.z)*SS_d.z;
                }
            }
        }

    return W.x + W.y + W.z;
}

double CoreVirial_Cicc(struct point r[], struct point rho[], int i){

    struct point W;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kQiQj, kQiChij, Aij, Cij, Dij;
    struct point CC_d, CS_d;
    double CC_r, CC_r2, CC_r3, CC_r8, CC_r10;
    double CS_r, CS_r3;

    CS_d.x = (r[i].x - rho[i].x);
    CS_d.y = (r[i].y - rho[i].y);
    CS_d.z = (r[i].z - rho[i].z);

    //    CS_d = Distance(r[i], rho[i]);
    //    CS_d = d_rirhoj(r[i], rho[i], r[i]);

    if (DEBUG_FLAG && _D_DIST) {
        printf("r[%d] = (%.4e, %.4e, %.4e), rho[%d] = (%.4e, %.4e, %.4e), CS_d = (%.4e, %.4e, %.4e)\n", i, r[i].x, r[i].y, r[i].z, i, rho[i].x, rho[i].y, rho[i].z, CS_d.x, CS_d.y, CS_d.z);
    }

    W.x = (-ki*CS_d.x)*CS_d.x;
    W.y = (-ki*CS_d.y)*CS_d.x;
    W.z = (-ki*CS_d.z)*CS_d.x;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Core-Core interactions
            //            CC_d = d_rirj(r[i], r[j]);
            CC_d = Distance(r[i], r[j]);
            if (DEBUG_FLAG && _D_DIST) printf("CC_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, CC_d.x, CC_d.y, CC_d.z);

            CC_r = mod(CC_d);
            if (DEBUG_FLAG && _D_DIST) printf("CC_r[%d][%d] = %.4e\n", i, j, CC_r);

            if (CC_r<LR_RCUT || R_CUT == -1) {

                kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];

                CC_r2 = CC_r*CC_r;
                CC_r3 = CC_r*CC_r2;

                W.x += (kQiQj/CC_r3*CC_d.x + kQiQj*LR_DVCUT/CC_r*CC_d.x + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.x)*CC_d.x;
                W.y += (kQiQj/CC_r3*CC_d.y + kQiQj*LR_DVCUT/CC_r*CC_d.y + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.y)*CC_d.y;
                W.z += (kQiQj/CC_r3*CC_d.z + kQiQj*LR_DVCUT/CC_r*CC_d.z + kQiQj*LR_DDVCUT*(CC_r - LR_RCUT)/CC_r*CC_d.z)*CC_d.z;
            }

            if (CC_r<RCUT || R_CUT == -1) {

                CC_r2 = CC_r*CC_r;

                Aij = A[indx_int];
                Cij = C[indx_int];
                Dij = D[indx_int];

                CC_r8 = CC_r2*CC_r2*CC_r2*CC_r2;
                CC_r10 = CC_r8*CC_r2;

                W.x += (Aij/LAMBDA*CC_d.x/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.x - 8.*Dij/CC_r10*CC_d.x + C_DVCUT[indx_int]/CC_r*CC_d.x + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.x)*CC_d.x;
                W.y += (Aij/LAMBDA*CC_d.y/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.y - 8.*Dij/CC_r10*CC_d.y + C_DVCUT[indx_int]/CC_r*CC_d.y + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.y)*CC_d.y;
                W.z += (Aij/LAMBDA*CC_d.z/CC_r*exp(-CC_r/LAMBDA) - 6.*Cij/CC_r8*CC_d.z - 8.*Dij/CC_r10*CC_d.z + C_DVCUT[indx_int]/CC_r*CC_d.z + C_DDVCUT[indx_int]*(CC_r - RCUT)/CC_r*CC_d.z)*CC_d.z;
            }

            //Core-Shell interactions
            //            CS_d = d_rirhoj(r[i], rho[j], r[j]);
            CS_d = Distance(r[i], rho[j]);
            if (DEBUG_FLAG && _D_DIST) printf("CS_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, CS_d.x, CS_d.y, CS_d.z);

            CS_r = mod(CS_d);
            if (DEBUG_FLAG && _D_DIST) printf("CS_r[%d][%d] = %.4e\n", i, j, CS_r);

            if (CS_r<LR_RCUT || R_CUT == -1) {

                kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                CS_r3 = CS_r*CS_r*CS_r;

                W.x += (kQiChij/CS_r3*CS_d.x + kQiChij*LR_DVCUT/CS_r*CS_d.x + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.x)*CS_d.x;
                W.y += (kQiChij/CS_r3*CS_d.y + kQiChij*LR_DVCUT/CS_r*CS_d.y + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.y)*CS_d.y;
                W.z += (kQiChij/CS_r3*CS_d.z + kQiChij*LR_DVCUT/CS_r*CS_d.z + kQiChij*LR_DDVCUT*(CS_r - LR_RCUT)/CS_r*CS_d.z)*CS_d.z;
            }
        }
    }

    return W.x + W.y + W.z;
}

double ShellVirial_Cicc(struct point rho[], struct point r[], int i){

    struct point W;

    int  j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj;
    struct point SS_d, SC_d;
    double SS_r, SS_r3;
    double SC_r, SC_r3;

    SC_d.x = (rho[i].x - r[i].x);
    SC_d.y = (rho[i].y - r[i].y);
    SC_d.z = (rho[i].z - r[i].z);

    //    SC_d = Distance(rho[i], r[i]);
    //    SC_d = d_rhoirj(rho[i], r[i], r[i]);

    if (DEBUG_FLAG && _D_DIST){

        printf("rho[%d] = (%.4e, %.4e, %.4e), r[%d] = (%.4e, %.4e, %.4e), SC_d = (%.4e, %.4e, %.4e)\n", i, rho[i].x, rho[i].y, rho[i].z, i, r[i].x, r[i].y, r[i].z, SC_d.x, SC_d.y, SC_d.z);
    }

    W.x = (-ki*SC_d.x)*SC_d.x;
    W.y = (-ki*SC_d.y)*SC_d.y;
    W.z = (-ki*SC_d.z)*SC_d.z;

    for(j=0; j<NPART; j++){

        if (i!=j){

            indx_j = INDX[j];
            indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

            //Shell-Shell interactions
            //            SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
            SS_d = Distance(rho[i], rho[j]);
            if (DEBUG_FLAG && _D_DIST) printf("SS_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, SS_d.x, SS_d.y, SS_d.z);

            SS_r = mod(SS_d);
            if (DEBUG_FLAG && _D_DIST) printf("SS_r[%d][%d] = %.4e\n", i, j, SS_r);

            if (SS_r<LR_RCUT || R_CUT == -1) {

                kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];

                SS_r3 = SS_r*SS_r*SS_r;

                W.x += (kChiiChij/SS_r3*SS_d.x + kChiiChij*LR_DVCUT/SS_r*SS_d.x + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.x)*SS_d.x;
                W.y += (kChiiChij/SS_r3*SS_d.y + kChiiChij*LR_DVCUT/SS_r*SS_d.y + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.y)*SS_d.y;
                W.z += (kChiiChij/SS_r3*SS_d.z + kChiiChij*LR_DVCUT/SS_r*SS_d.z + kChiiChij*LR_DDVCUT*(SS_r - LR_RCUT)/SS_r*SS_d.z)*SS_d.z;
            }

            //Shell-Core interactions
            //            SC_d = d_rhoirj(rho[i], r[i], r[j]);
            SC_d = Distance(rho[i], r[j]);
            if (DEBUG_FLAG && _D_DIST) printf("SC_d[%d][%d] = (%.4e, %.4e, %.4e)\n", i, j, SC_d.x, SC_d.y, SC_d.z);

            SC_r = mod(SC_d);
            if (DEBUG_FLAG && _D_DIST) printf("SC_r[%d][%d] = %.4e\n", i, j, SC_r);

            if (SC_r<LR_RCUT || R_CUT == -1) {

                kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                SC_r3 = SC_r*SC_r*SC_r;

                W.x += (kChiiQj/SC_r3*SC_d.x + kChiiQj*LR_DVCUT/SC_r*SC_d.x + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.x)*SC_d.x;
                W.y += (kChiiQj/SC_r3*SC_d.y + kChiiQj*LR_DVCUT/SC_r*SC_d.y + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.y)*SC_d.y;
                W.z += (kChiiQj/SC_r3*SC_d.z + kChiiQj*LR_DVCUT/SC_r*SC_d.z + kChiiQj*LR_DDVCUT*(SC_r - LR_RCUT)/SC_r*SC_d.z)*SC_d.z;
            }
        }
    }

    return W.x + W.y + W.z;
}

void Forces(struct point r[], struct point rho[]){

    int i, j;
    int indx_i, indx_j, indx_int, indx_g;

    struct point SC_d_ii, CC_d_ij, CS_d_ij, SC_d_ij, SS_d_ij;
    double CC_r, CC_r2, CC_r6, CC_r8, CS_r, CS_r2, SC_r, SC_r2, SS_r, SS_r2;

    double kQiQj, kQiChij, kChiiQj, kChiiChij, Aij, Cij, Dij;

    double W_ij;

    struct point CF_ij, SF_ij, CF_tot, SF_tot;

    EPOT = 0;
    PRESS = 0;
    NGCOUNT++;

    CF_tot.x = CF_tot.y = CF_tot.z = SF_tot.x = SF_tot.y = SF_tot.z = 0;

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

        SC_d_ii.x = rho[i].x - r[i].x;
        SC_d_ii.y = rho[i].y - r[i].y;
        SC_d_ii.z = rho[i].z - r[i].z;

        EPOT += K[indx_i]*(SC_d_ii.x*SC_d_ii.x + SC_d_ii.y*SC_d_ii.y + SC_d_ii.z*SC_d_ii.z);

        CF[i].x = K[indx_i]*SC_d_ii.x;
        CF[i].x = K[indx_i]*SC_d_ii.y;
        CF[i].x = K[indx_i]*SC_d_ii.z;

        SF[i].x = -CF[i].x;
        SF[i].x = -CF[i].x;
        SF[i].x = -CF[i].x;
    }

    for (i=0; i<(NPART-1); i++) {

//        printf("i = %d", i);

        indx_i = INDX[i];

        for (j=i+1; j<NPART; j++) {

//            printf("\t\tj = %d\n", j);

            indx_j = INDX[j];
            indx_int = indx_i + indx_j;

            CC_d_ij.x = r[i].x - r[j].x;
            CC_d_ij.x -= LBOX*nearbyint(CC_d_ij.x/LBOX);
            CC_d_ij.y = r[i].y - r[j].y;
            CC_d_ij.y -= LBOX*nearbyint(CC_d_ij.y/LBOX);
            CC_d_ij.z = r[i].z - r[j].z;
            CC_d_ij.z -= LBOX*nearbyint(CC_d_ij.z/LBOX);

            CC_r = sqrt(CC_d_ij.x*CC_d_ij.x + CC_d_ij.y*CC_d_ij.y + CC_d_ij.z*CC_d_ij.z);

            if (CC_r < RCUT) {

                CC_r2 = CC_r*CC_r;
                CC_r6 = CC_r2*CC_r2*CC_r2;
                CC_r8 = CC_r6*CC_r2;

                kQiQj = _COULOMB*Q[indx_i]*Q[indx_j];
                Cij = C[indx_int];
                Dij = D[indx_int];

                EPOT += kQiQj*(erfc(ALPHA*CC_r)/CC_r - EW_VCUT) - Cij/CC_r6 - Dij/CC_r8 - C_VCUT[indx_int] - (kQiQj*EW_DVCUT + C_DVCUT[indx_int])*(CC_r - RCUT) - .5*(kQiQj*EW_DDVCUT + C_DDVCUT[indx_int])*(CC_r - RCUT)*(CC_r - RCUT);

                W_ij = -kQiQj*(M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*CC_r2) + erfc(ALPHA*CC_r)/CC_r) + 6.*Cij/CC_r6 + 8.*Dij/CC_r8 - (kQiQj*EW_DVCUT + C_DVCUT[indx_int])*CC_r - (kQiQj*EW_DDVCUT + C_DDVCUT[indx_int])*(CC_r - RCUT)*CC_r;

                PRESS += W_ij;

                CF_ij.x = -W_ij/CC_r2*CC_d_ij.x;
                CF_ij.y = -W_ij/CC_r2*CC_d_ij.y;
                CF_ij.z = -W_ij/CC_r2*CC_d_ij.z;

                CF[i].x += CF_ij.x;
                CF[i].y += CF_ij.y;
                CF[i].z += CF_ij.z;

                CF[j].x -= CF_ij.x;
                CF[j].y -= CF_ij.y;
                CF[j].z -= CF_ij.z;
            }

            if (CC_r < .5*LBOX) {

                indx_g = (int) floor(CC_r/DG);
                GCOUNT[indx_int][indx_g] += 2;
            }

            CS_d_ij.x = r[i].x - rho[j].x;
            CS_d_ij.x -= LBOX*nearbyint(CS_d_ij.x/LBOX);
            CS_d_ij.y = r[i].y - rho[j].y;
            CS_d_ij.y -= LBOX*nearbyint(CS_d_ij.y/LBOX);
            CS_d_ij.z = r[i].z - rho[j].z;
            CS_d_ij.z -= LBOX*nearbyint(CS_d_ij.z/LBOX);

            CS_r = sqrt(CS_d_ij.x*CS_d_ij.x + CS_d_ij.y*CS_d_ij.y + CS_d_ij.z*CS_d_ij.z);

            if (CS_r < RCUT) {

                CS_r2 = CS_r*CS_r;

                kQiChij = _COULOMB*Q[indx_i]*CHI[indx_j];

                EPOT += kQiChij*(erfc(ALPHA*CS_r)/CS_r - EW_VCUT - EW_DVCUT*(CS_r - RCUT) - .5*EW_DDVCUT*(CS_r - RCUT)*(CS_r - RCUT));

                W_ij = -kQiChij*(M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*CS_r2) + erfc(ALPHA*CS_r)/CS_r + EW_DVCUT*CS_r + EW_DDVCUT*(CS_r - RCUT)*CS_r);

                PRESS += W_ij;

                CF_ij.x = -W_ij/CS_r2*CS_d_ij.x;
                CF_ij.y = -W_ij/CS_r2*CS_d_ij.y;
                CF_ij.z = -W_ij/CS_r2*CS_d_ij.z;

                CF[i].x += CF_ij.x;
                CF[i].y += CF_ij.y;
                CF[i].z += CF_ij.z;

                CF[j].x -= CF_ij.x;
                CF[j].y -= CF_ij.y;
                CF[j].z -= CF_ij.z;
            }

            SC_d_ij.x = rho[i].x - r[j].x;
            SC_d_ij.x -= LBOX*nearbyint(SC_d_ij.x/LBOX);
            SC_d_ij.y = rho[i].y - r[j].y;
            SC_d_ij.y -= LBOX*nearbyint(SC_d_ij.y/LBOX);
            SC_d_ij.z = rho[i].z - r[j].z;
            SC_d_ij.z -= LBOX*nearbyint(SC_d_ij.z/LBOX);

            SC_r = sqrt(SC_d_ij.x*SC_d_ij.x + SC_d_ij.y*SC_d_ij.y + SC_d_ij.z*SC_d_ij.z);

            if (SC_r < RCUT) {

                SC_r2 = SC_r*SC_r;

                kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];

                EPOT += kChiiQj*(erfc(ALPHA*SC_r)/SC_r - EW_VCUT - EW_DVCUT*(SC_r - RCUT) - .5*EW_DDVCUT*(SC_r - RCUT)*(SC_r - RCUT));

                W_ij = -kChiiQj*(M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*SC_r2) + erfc(ALPHA*SC_r)/SC_r + EW_DVCUT*SC_r + EW_DDVCUT*(SC_r - RCUT)*SC_r);

                PRESS += W_ij;

                SF_ij.x = -W_ij/SC_r2*SC_d_ij.x;
                SF_ij.y = -W_ij/SC_r2*SC_d_ij.y;
                SF_ij.z = -W_ij/SC_r2*SC_d_ij.z;

                SF[i].x += SF_ij.x;
                SF[i].y += SF_ij.y;
                SF[i].z += SF_ij.z;

                SF[j].x -= SF_ij.x;
                SF[j].y -= SF_ij.y;
                SF[j].z -= SF_ij.z;
            }

            SS_d_ij.x = rho[i].x - rho[j].x;
            SS_d_ij.x -= LBOX*nearbyint(SS_d_ij.x/LBOX);
            SS_d_ij.y = rho[i].y - rho[j].y;
            SS_d_ij.y -= LBOX*nearbyint(SS_d_ij.y/LBOX);
            SS_d_ij.z = rho[i].z - rho[j].z;
            SS_d_ij.z -= LBOX*nearbyint(SS_d_ij.z/LBOX);

            SS_r = sqrt(SS_d_ij.x*SS_d_ij.x + SS_d_ij.y*SS_d_ij.y + SS_d_ij.z*SS_d_ij.z);

            if (SS_r < RCUT) {

                SS_r2 = SS_r*SS_r;

                kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];
                Aij = A[indx_int];

                EPOT += kChiiChij*(erfc(ALPHA*SS_r)/SS_r - EW_VCUT) + Aij*exp(-SS_r/LAMBDA) - S_VCUT[indx_int] - (kChiiChij*EW_DVCUT + S_DVCUT[indx_int])*(SS_r - RCUT) - .5*(kChiiChij*EW_DDVCUT + S_DDVCUT[indx_int])*(SS_r - RCUT)*(SS_r - RCUT);

                W_ij = -kChiiChij*(M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*SS_r2) + erfc(ALPHA*SS_r)/SS_r) - Aij/LAMBDA*SS_r*exp(-SS_r/LAMBDA) - (kChiiChij*EW_DVCUT + S_DVCUT[indx_int])*SS_r - (kChiiChij*EW_DDVCUT + S_DDVCUT[indx_int])*(SS_r - RCUT)*SS_r;

                PRESS += W_ij;

                SF_ij.x = -W_ij/SS_r2*SS_d_ij.x;
                SF_ij.y = -W_ij/SS_r2*SS_d_ij.y;
                SF_ij.z = -W_ij/SS_r2*SS_d_ij.z;

                SF[i].x += SF_ij.x;
                SF[i].y += SF_ij.y;
                SF[i].z += SF_ij.z;

                SF[j].x -= SF_ij.x;
                SF[j].y -= SF_ij.y;
                SF[j].z -= SF_ij.z;
            }
        }
    }

    for (i=0; i<NPART; i++) {

        printf("CF[%d] = (%lf, %lf, %lf)\tSF[%d] = (%lf, %lf, %lf)\n", i, CF[i].x, CF[i].y, CF[i].z, i, SF[i].x, SF[i].y, SF[i].z);

        CF_tot.x += CF[i].x;
        CF_tot.y += CF[i].y;
        CF_tot.z += CF[i].z;

        SF_tot.x += SF[i].x;
        SF_tot.y += SF[i].y;
        SF_tot.z += SF[i].z;
    }

//    printf("\n\nCF_tot = (%.4e, %.4e, %.4e)\tSf_tot = (%.4e, %.4e, %.4e)\n\n", CF_tot.x, CF_tot.y, CF_tot.z, SF_tot.x, SF_tot.y, SF_tot.z);
}
