//
//  ConstTensor.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "ConstTensor.h"

struct tensor ConstTens_Jac(struct point rho[], struct point r[], int i, int k){

    struct tensor W;

    int j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj, Aij;

    struct point SS_d, SC_d;
    double SS_r, SS_r2, SS_r3, m3kchichi_SSr5i, kchichi3_SSr5i, mSSr2i_expmSSr, SSr2i_expmSSr;
    double SC_r, SC_r2, SC_r3, m3kchiQ_SCr5i;

    if (i==k){

        W.fx.x = W.fy.y = W.fz.z = -ki;
        W.fx.y = W.fx.z = W.fy.x = W.fy.z = W.fz.x = W.fz.y = 0.;

        for (j=0; j<NPART; j++){

            if (i!=j){

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

//                SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (SS_r < LR_RCUT || R_CUT == -1){

                    SS_r2 = SS_r*SS_r;
                    SS_r3 = SS_r2*SS_r;

                    kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];
                    m3kchichi_SSr5i = -3.*kChiiChij/SS_r2/SS_r2/SS_r;

                    W.fx.x += m3kchichi_SSr5i*(SS_d.x*SS_d.x-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
                    W.fx.y += m3kchichi_SSr5i*(SS_d.x*SS_d.y) - kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.y + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.y;
                    W.fx.z += m3kchichi_SSr5i*(SS_d.x*SS_d.z) - kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.z + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.z;
                    W.fy.x += m3kchichi_SSr5i*(SS_d.y*SS_d.x) - kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.x + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.x;
                    W.fy.y += m3kchichi_SSr5i*(SS_d.y*SS_d.y-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
                    W.fy.z += m3kchichi_SSr5i*(SS_d.y*SS_d.z) - kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.z + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.z;
                    W.fz.x += m3kchichi_SSr5i*(SS_d.z*SS_d.x) - kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.x + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.x;
                    W.fz.y += m3kchichi_SSr5i*(SS_d.z*SS_d.y) - kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.y + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.y;
                    W.fz.z += m3kchichi_SSr5i*(SS_d.z*SS_d.z-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
                }

                if (SS_r < RCUT || R_CUT == -1){

                    SS_r2 = SS_r*SS_r;
                    SS_r3 = SS_r2*SS_r;

                    Aij = A[indx_int];
                    mSSr2i_expmSSr = -Aij/LAMBDA/SS_r2*exp(-SS_r/LAMBDA);

                    W.fx.x += mSSr2i_expmSSr*((SS_d.x*SS_d.x-SS_r2)/SS_r + (SS_d.x*SS_d.x)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) + S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
                    W.fx.y += mSSr2i_expmSSr*((SS_d.x*SS_d.y)/SS_r + (SS_d.x*SS_d.y)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.x*SS_d.y + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.x*SS_d.y;
                    W.fx.z += mSSr2i_expmSSr*((SS_d.x*SS_d.z)/SS_r + (SS_d.x*SS_d.z)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.x*SS_d.z + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.x*SS_d.z;
                    W.fy.x += mSSr2i_expmSSr*((SS_d.y*SS_d.x)/SS_r + (SS_d.y*SS_d.x)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.y*SS_d.x + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.y*SS_d.x;
                    W.fy.y += mSSr2i_expmSSr*((SS_d.y*SS_d.y-SS_r2)/SS_r + (SS_d.y*SS_d.y)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) + S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
                    W.fy.z += mSSr2i_expmSSr*((SS_d.y*SS_d.z)/SS_r + (SS_d.y*SS_d.z)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.y*SS_d.z + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.y*SS_d.z;
                    W.fz.x += mSSr2i_expmSSr*((SS_d.z*SS_d.x)/SS_r + (SS_d.z*SS_d.x)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.z*SS_d.x + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.z*SS_d.x;
                    W.fz.y += mSSr2i_expmSSr*((SS_d.z*SS_d.y)/SS_r + (SS_d.z*SS_d.y)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*SS_d.z*SS_d.y + S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.z*SS_d.y;
                    W.fz.z += mSSr2i_expmSSr*((SS_d.z*SS_d.z-SS_r2)/SS_r + (SS_d.z*SS_d.z)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) + S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
                }

//                SC_d = d_rhoirj(rho[i], r[i], r[j]);
                SC_d = Distance(rho[i], r[j]);
                SC_r = mod(SC_d);

                if (SC_r < LR_RCUT || R_CUT == -1){

                    SC_r2 = SC_r*SC_r;
                    SC_r3 = SC_r2*SC_r;

                    kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];
                    m3kchiQ_SCr5i = -3.*kChiiQj/SC_r2/SC_r2/SC_r;

                    W.fx.x += m3kchiQ_SCr5i*(SC_d.x*SC_d.x - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.y*SC_d.y + SC_d.z*SC_d.z) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.y*SC_d.y + SC_d.z*SC_d.z));
                    W.fx.y += m3kchiQ_SCr5i*SC_d.x*SC_d.y - kChiiQj*LR_DVCUT/SC_r3*SC_d.x*SC_d.y + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.x*SC_d.y;
                    W.fx.z += m3kchiQ_SCr5i*SC_d.x*SC_d.z - kChiiQj*LR_DVCUT/SC_r3*SC_d.x*SC_d.z + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.x*SC_d.z;
                    W.fy.x += m3kchiQ_SCr5i*SC_d.y*SC_d.x - kChiiQj*LR_DVCUT/SC_r3*SC_d.y*SC_d.x + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.y*SC_d.x;
                    W.fy.y += m3kchiQ_SCr5i*(SC_d.y*SC_d.y - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.z*SC_d.z) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.z*SC_d.z));
                    W.fy.z += m3kchiQ_SCr5i*SC_d.y*SC_d.z - kChiiQj*LR_DVCUT/SC_r3*SC_d.y*SC_d.z + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.y*SC_d.z;
                    W.fz.x += m3kchiQ_SCr5i*SC_d.z*SC_d.x - kChiiQj*LR_DVCUT/SC_r3*SC_d.z*SC_d.x + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.z*SC_d.x;
                    W.fz.y += m3kchiQ_SCr5i*SC_d.z*SC_d.y - kChiiQj*LR_DVCUT/SC_r3*SC_d.z*SC_d.y + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.z*SC_d.y;
                    W.fz.z += m3kchiQ_SCr5i*(SC_d.z*SC_d.z - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.y*SC_d.y) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.y*SC_d.y));
                }
            }
        }

    }else{

        W.fx.x = W.fx.y = W.fx.z = W.fy.x = W.fy.y = W.fy.z = W.fz.x = W.fz.y = W.fz.z = 0.;

        j=k;

        indx_j = INDX[j];
        indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

//        SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
        SS_d = Distance(rho[i], rho[j]);
        SS_r = mod(SS_d);

        if (SS_r < LR_RCUT || R_CUT == -1){

            SS_r2 = SS_r*SS_r;
            SS_r3 = SS_r2*SS_r;

            kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];
            kchichi3_SSr5i = kChiiChij*3./SS_r2/SS_r2/SS_r;

            W.fx.x = kchichi3_SSr5i*(SS_d.x*SS_d.x-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
            W.fx.y = kchichi3_SSr5i*SS_d.x*SS_d.y + kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.y - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.y;
            W.fx.z = kchichi3_SSr5i*SS_d.x*SS_d.z + kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.z - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.z;
            W.fy.x = kchichi3_SSr5i*SS_d.y*SS_d.x + kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.x - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.x;
            W.fy.y = kchichi3_SSr5i*(SS_d.y*SS_d.y-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
            W.fy.z = kchichi3_SSr5i*SS_d.y*SS_d.z + kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.z - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.z;
            W.fz.x = kchichi3_SSr5i*SS_d.z*SS_d.x + kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.x - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.x;
            W.fz.y = kchichi3_SSr5i*SS_d.z*SS_d.y + kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.y - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.y;
            W.fz.z = kchichi3_SSr5i*(SS_d.z*SS_d.z-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
        }

        if (SS_r < RCUT || R_CUT == -1){

            SS_r2 = SS_r*SS_r;
            SS_r3 = SS_r2*SS_r;

            Aij = A[indx_int];
            SSr2i_expmSSr = Aij/LAMBDA/SS_r2*exp(-SS_r/LAMBDA);

            W.fx.x = SSr2i_expmSSr*((SS_d.x*SS_d.x-SS_r2)/SS_r + (SS_d.x*SS_d.x)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) - S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
            W.fx.y = SSr2i_expmSSr*((SS_d.x*SS_d.y)/SS_r + (SS_d.x*SS_d.y)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.x*SS_d.y - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.x*SS_d.y;
            W.fx.z = SSr2i_expmSSr*((SS_d.x*SS_d.z)/SS_r + (SS_d.x*SS_d.z)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.x*SS_d.z - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.x*SS_d.z;
            W.fy.x = SSr2i_expmSSr*((SS_d.y*SS_d.x)/SS_r + (SS_d.y*SS_d.x)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.y*SS_d.x - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.y*SS_d.x;
            W.fy.y = SSr2i_expmSSr*((SS_d.y*SS_d.y-SS_r2)/SS_r + (SS_d.y*SS_d.y)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) - S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
            W.fy.z = SSr2i_expmSSr*((SS_d.y*SS_d.z)/SS_r + (SS_d.y*SS_d.z)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.y*SS_d.z - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.y*SS_d.z;
            W.fz.x = SSr2i_expmSSr*((SS_d.z*SS_d.x)/SS_r + (SS_d.z*SS_d.x)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.z*SS_d.x - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.z*SS_d.x;
            W.fz.y = SSr2i_expmSSr*((SS_d.z*SS_d.y)/SS_r + (SS_d.z*SS_d.y)/LAMBDA) + S_DVCUT[indx_int]/SS_r3*SS_d.z*SS_d.y - S_DDVCUT[indx_int]*RCUT/SS_r3*SS_d.z*SS_d.y;
            W.fz.z = SSr2i_expmSSr*((SS_d.z*SS_d.z-SS_r2)/SS_r + (SS_d.z*SS_d.z)/LAMBDA) - S_DVCUT[indx_int]/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) - S_DDVCUT[indx_int]*(1-RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
        }
    }
    return W;
}

struct tensor ConstTens_Cicc(struct point rho[], struct point r[], int i, int k){

    struct tensor W;

    int j, indx_i = INDX[i], indx_j, indx_int;
    double ki = K[indx_i], kChiiChij, kChiiQj;

    struct point SS_d, SC_d;
    double SS_r, SS_r2, SS_r3, m3kchichi_SSr5i, kchichi3_SSr5i;
    double SC_r, SC_r2, SC_r3, m3kchiQ_SCr5i;

    if (i==k){

        W.fx.x = W.fy.y = W.fz.z = -ki;
        W.fx.y = W.fx.z = W.fy.x = W.fy.z = W.fz.x = W.fz.y = 0.;

        for (j=0; j<NPART; j++){

            if (i!=j){

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

//                SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (SS_r <LR_RCUT || R_CUT == -1){

                    SS_r2 = SS_r*SS_r;
                    SS_r3 = SS_r2*SS_r;

                    kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];
                    m3kchichi_SSr5i = -3.*kChiiChij/SS_r2/SS_r2/SS_r;

                    W.fx.x += m3kchichi_SSr5i*(SS_d.x*SS_d.x-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
                    W.fx.y += m3kchichi_SSr5i*(SS_d.x*SS_d.y) - kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.y + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.y;
                    W.fx.z += m3kchichi_SSr5i*(SS_d.x*SS_d.z) - kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.z + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.z;
                    W.fy.x += m3kchichi_SSr5i*(SS_d.y*SS_d.x) - kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.x + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.x;
                    W.fy.y += m3kchichi_SSr5i*(SS_d.y*SS_d.y-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
                    W.fy.z += m3kchichi_SSr5i*(SS_d.y*SS_d.z) - kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.z + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.z;
                    W.fz.x += m3kchichi_SSr5i*(SS_d.z*SS_d.x) - kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.x + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.x;
                    W.fz.y += m3kchichi_SSr5i*(SS_d.z*SS_d.y) - kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.y + kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.y;
                    W.fz.z += m3kchichi_SSr5i*(SS_d.z*SS_d.z-SS_r2/3.) + kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) + kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
                }

//                SC_d = d_rhoirj(rho[i], r[i], r[j]);
                SC_d = Distance(rho[i], r[j]);
                SC_r = mod(SC_d);

                if (SC_r <LR_RCUT || R_CUT == -1){

                    SC_r2 = SC_r*SC_r;
                    SC_r3 = SC_r2*SC_r;

                    kChiiQj = _COULOMB*CHI[indx_i]*Q[indx_j];
                    m3kchiQ_SCr5i = -3.*kChiiQj/SC_r2/SC_r2/SC_r;

                    W.fx.x += m3kchiQ_SCr5i*(SC_d.x*SC_d.x - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.y*SC_d.y + SC_d.z*SC_d.z) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.y*SC_d.y + SC_d.z*SC_d.z));
                    W.fx.y += m3kchiQ_SCr5i*(SC_d.x*SC_d.y) - kChiiQj*LR_DVCUT/SC_r3*SC_d.x*SC_d.y + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.x*SC_d.y;
                    W.fx.z += m3kchiQ_SCr5i*(SC_d.x*SC_d.z) - kChiiQj*LR_DVCUT/SC_r3*SC_d.x*SC_d.z + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.x*SC_d.z;
                    W.fy.x += m3kchiQ_SCr5i*(SC_d.y*SC_d.x) - kChiiQj*LR_DVCUT/SC_r3*SC_d.y*SC_d.x + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.y*SC_d.x;
                    W.fy.y += m3kchiQ_SCr5i*(SC_d.y*SC_d.y - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.z*SC_d.z) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.z*SC_d.z));
                    W.fy.z += m3kchiQ_SCr5i*(SC_d.y*SC_d.z) - kChiiQj*LR_DVCUT/SC_r3*SC_d.y*SC_d.z + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.y*SC_d.z;
                    W.fz.x += m3kchiQ_SCr5i*(SC_d.z*SC_d.x) - kChiiQj*LR_DVCUT/SC_r3*SC_d.z*SC_d.x + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.z*SC_d.x;
                    W.fz.y += m3kchiQ_SCr5i*(SC_d.z*SC_d.y) - kChiiQj*LR_DVCUT/SC_r3*SC_d.z*SC_d.y + kChiiQj*LR_DDVCUT*LR_RCUT/SC_r3*SC_d.z*SC_d.y;
                    W.fz.z += m3kchiQ_SCr5i*(SC_d.z*SC_d.z - SC_r2/3.) + kChiiQj*LR_DVCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.y*SC_d.y) + kChiiQj*LR_DDVCUT*(1-LR_RCUT/SC_r3*(SC_d.x*SC_d.x + SC_d.y*SC_d.y));
                }
            }
        }

    }else{

        W.fx.x = W.fx.y = W.fx.z = W.fy.x = W.fy.y = W.fy.z = W.fz.x = W.fz.y = W.fz.z = 0.;

        j=k;

        indx_j = INDX[j];
        indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

//        SS_d = d_rhoirhoj(rho[i], r[i], rho[j], r[j]);
        SS_d = Distance(rho[i], rho[j]);
        SS_r = mod(SS_d);

        if (SS_r < LR_RCUT || R_CUT == -1){

            SS_r2 = SS_r*SS_r;
            SS_r3 = SS_r2*SS_r;

            kChiiChij = _COULOMB*CHI[indx_i]*CHI[indx_j];
            kchichi3_SSr5i = kChiiChij*3./SS_r2/SS_r2/SS_r;

            W.fx.x = kchichi3_SSr5i*(SS_d.x*SS_d.x-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.y*SS_d.y + SS_d.z*SS_d.z));
            W.fx.y = kchichi3_SSr5i*SS_d.x*SS_d.y + kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.y - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.y;
            W.fx.z = kchichi3_SSr5i*SS_d.x*SS_d.z + kChiiChij*LR_DVCUT/SS_r3*SS_d.x*SS_d.z - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.x*SS_d.z;
            W.fy.x = kchichi3_SSr5i*SS_d.y*SS_d.x + kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.x - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.x;
            W.fy.y = kchichi3_SSr5i*(SS_d.y*SS_d.y-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.z*SS_d.z));
            W.fy.z = kchichi3_SSr5i*SS_d.y*SS_d.z + kChiiChij*LR_DVCUT/SS_r3*SS_d.y*SS_d.z - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.y*SS_d.z;
            W.fz.x = kchichi3_SSr5i*SS_d.z*SS_d.x + kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.x - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.x;
            W.fz.y = kchichi3_SSr5i*SS_d.z*SS_d.y + kChiiChij*LR_DVCUT/SS_r3*SS_d.z*SS_d.y - kChiiChij*LR_DDVCUT*LR_RCUT/SS_r3*SS_d.z*SS_d.y;
            W.fz.z = kchichi3_SSr5i*(SS_d.z*SS_d.z-SS_r2/3.) - kChiiChij*LR_DVCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y) - kChiiChij*LR_DDVCUT*(1-LR_RCUT/SS_r3*(SS_d.x*SS_d.x + SS_d.y*SS_d.y));
        }
    }

    return W;
}

struct tensor ConstTens_WCA(struct point rho[], struct point r[], int k, int i) {

//Here k is the index of the constraint (one component of the force) and i is the index in the derivative
    struct tensor W;

    struct point SS_d;
    double SS_r;
    double rCUT, lround, Rround, R2round;
    int j, indx_i = INDX[i], indx_j, indx_int;
    struct point CC_d;
    double CC_r, rinv, r2inv, r6inv, rc2inv, rc6inv, ljatrc, forcelj, fpair, ljrc;
    double LJ_sigma6, LJ_sigma12;
    double sigma_r6, sigma_r12, dU, dU2;


    if (i==k) {

        W.fx.x = W.fx.y = W.fx.z = 0.;
        W.fy.x = W.fy.y = W.fy.z = 0.;
        W.fz.x = W.fz.y = W.fz.z = 0.;

        //for (j=0; j<NPART; j++){
        int p;
        int neighlist[1000];
        List_Of_Neighs(k,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {

            j = neighlist[p];

            if (j!=k) {

                indx_j = INDX[j];
                indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

                rCUT = LJRCUT[indx_i][indx_j];
                lround = LJLROUND[indx_i][indx_j];
                LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
                LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);
                SS_d = Distance(rho[i], rho[j]);
                SS_r = mod(SS_d);

                if (SS_r <= rCUT) {



                    rinv = 1.0/SS_r;
                    r2inv = rinv*rinv;
                    r6inv = r2inv*r2inv*r2inv;
                    rc2inv = 1.0/(rCUT*rCUT);
                    rc6inv = rc2inv*rc2inv*rc2inv;
                    sigma_r6 = LJ_sigma6*r6inv;
                    sigma_r12 = LJ_sigma12*r6inv*r6inv;
                    ljatrc = rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
                    //ljatrc = 4. * LJEPS[indx_i][indx_j] * rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
                    dU = 24.0 * LJEPS[indx_i][indx_j]*rinv * (-2.0*sigma_r12 + sigma_r6);
                    dU2 = 24.0 * LJEPS[indx_i][indx_j]*r2inv * (26.0*sigma_r12 - 7.0*sigma_r6);

                    //rounding the force

                    if (SS_r > (rCUT-lround)) {
                        Rround = (SS_r - rCUT + lround)/lround;
                        //forcelj = 24.0
                        //          * (Rround * (Rround-1.0)/lround * (LJEPS[indx_i][indx_j] * r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + 0.25*ljatrc) //0.25 missing?
                        //          + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * LJEPS[indx_i][indx_j] * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
                        dU = 24.0 * LJEPS[indx_i][indx_j] * (Rround * (Rround-1.0)/lround * (sigma_r12 - sigma_r6 - ljatrc)
                             + (1.0 + Rround * Rround * (2.0*Rround - 3.0))*rinv * (-2.0*sigma_r12 + sigma_r6));

                        dU2 = 24.0 * LJEPS[indx_i][indx_j] * ((2.0*Rround - 1.0)/(lround*lround) * (sigma_r12 - sigma_r6 - ljatrc)
                              + 12.0*(Rround * (Rround-1.0))/lround*rinv * (-2.0*sigma_r12 + sigma_r6)
                                + (1.0 + Rround * Rround * (2.0*Rround - 3.0))*r2inv * (26.0*sigma_r12 - 7.0*sigma_r6));

                        //forcelj = 24.0 * LJEPS[indx_i][indx_j]
                                  // * (Rround * (Rround-1.0)/lround * (r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + ljatrc)
                                  //    + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
                    }



                    W.fx.x += (- SS_d.x*SS_d.x*r2inv*dU2 - (SS_r - SS_d.x*SS_d.x*rinv)*r2inv*dU);
                    W.fx.y += (- SS_d.x*SS_d.y*r2inv*dU2 + (SS_d.x*SS_d.y*rinv)*r2inv*dU);
                    W.fx.z += (- SS_d.x*SS_d.z*r2inv*dU2 + (SS_d.x*SS_d.z*rinv)*r2inv*dU);
                    W.fy.x += (- SS_d.y*SS_d.x*r2inv*dU2 + (SS_d.y*SS_d.x*rinv)*r2inv*dU);
                    W.fy.y += (- SS_d.y*SS_d.y*r2inv*dU2 - (SS_r - SS_d.y*SS_d.y*rinv*r2inv)*dU);
                    W.fy.z += (- SS_d.y*SS_d.z*r2inv*dU2 + (SS_d.y*SS_d.z*rinv)*r2inv*dU);
                    W.fz.x += (- SS_d.z*SS_d.x*r2inv*dU2 + (SS_d.z*SS_d.x*rinv)*r2inv*dU);
                    W.fz.y += (- SS_d.z*SS_d.y*r2inv*dU2 + (SS_d.z*SS_d.y*rinv)*r2inv*dU);
                    W.fz.z += (- SS_d.z*SS_d.z*r2inv*dU2 - (SS_r - SS_d.z*SS_d.z*rinv)*r2inv*dU);
                }
            }
        }

    } else {

        W.fx.x = W.fx.y = W.fx.z = W.fy.x = W.fy.y = W.fy.z = W.fz.x = W.fz.y = W.fz.z = 0.;

        j=k;


        indx_j = INDX[j];
        indx_int = indx_i+indx_j; //indx_int = 0 -> ANAN, indx_int = 1 -> ANACAT, indx_int = 2 -> CATCAT

        rCUT = LJRCUT[indx_i][indx_j];
        lround = LJLROUND[indx_i][indx_j];
        LJ_sigma6 = pow((double)LJSIGMA[indx_i][indx_j],6.);
        LJ_sigma12 = pow((double)LJSIGMA[indx_i][indx_j],12.);
        SS_d = Distance(rho[i], rho[j]);
        SS_r = mod(SS_d);

        if (SS_r <= rCUT) {



            rinv = 1.0/SS_r;
            r2inv = rinv*rinv;
            r6inv = r2inv*r2inv*r2inv;
            rc2inv = 1.0/(rCUT*rCUT);
            rc6inv = rc2inv*rc2inv*rc2inv;
            sigma_r6 = LJ_sigma6*r6inv;
            sigma_r12 = LJ_sigma12*r6inv*r6inv;
            ljatrc = rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
            //ljatrc = 4. * LJEPS[indx_i][indx_j] * rc6inv * (LJ_sigma12*rc6inv - LJ_sigma6); // computing the shift
            dU = 24.0 * LJEPS[indx_i][indx_j]*rinv * (-2.0*sigma_r12 + sigma_r6);
            dU2 = 24.0 * LJEPS[indx_i][indx_j]*r2inv * (26.0*sigma_r12 - 7.0*sigma_r6);

            //rounding the force

            if (SS_r > (rCUT-lround)) {
                Rround = (SS_r - rCUT + lround)/lround;
                //forcelj = 24.0
                //          * (Rround * (Rround-1.0)/lround * (LJEPS[indx_i][indx_j] * r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + 0.25*ljatrc) //0.25 missing?
                //          + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * LJEPS[indx_i][indx_j] * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
                dU = 24.0 * LJEPS[indx_i][indx_j] * (Rround * (Rround-1.0)/lround * (sigma_r12 - sigma_r6 - ljatrc)
                     + (1.0 + Rround * Rround * (2.0*Rround - 3.0))*rinv * (-2.0*sigma_r12 + sigma_r6));

                dU2 = 24.0 * LJEPS[indx_i][indx_j] * ((2.0*Rround - 1.0)/(lround*lround) * (sigma_r12 - sigma_r6 - ljatrc)
                      + 12.0*(Rround * (Rround-1.0))/lround*rinv * (-2.0*sigma_r12 + sigma_r6)
                        + (1.0 + Rround * Rround * (2.0*Rround - 3.0))*r2inv * (26.0*sigma_r12 - 7.0*sigma_r6));

                //forcelj = 24.0 * LJEPS[indx_i][indx_j]
                          // * (Rround * (Rround-1.0)/lround * (r6inv * (LJ_sigma6 - LJ_sigma12*r6inv) + ljatrc)
                          //    + (1.0 + Rround * Rround * (2.0*Rround - 3.0)) * r6inv/CC_r * (2.0*LJ_sigma12*r6inv - LJ_sigma6));
            }

            W.fx.x = SS_d.x*SS_d.x*r2inv*dU2 + (SS_r - SS_d.x*SS_d.x*rinv)*r2inv*dU;
            W.fx.y = SS_d.x*SS_d.y*r2inv*dU2 - (SS_d.x*SS_d.y*rinv)*r2inv*dU;
            W.fx.z = SS_d.x*SS_d.z*r2inv*dU2 - (SS_d.x*SS_d.z*rinv)*r2inv*dU;
            W.fy.x = SS_d.y*SS_d.x*r2inv*dU2 - (SS_d.y*SS_d.x*rinv)*r2inv*dU;
            W.fy.y = SS_d.y*SS_d.y*r2inv*dU2 + (SS_r - SS_d.y*SS_d.y*rinv*r2inv)*dU;
            W.fy.z = SS_d.y*SS_d.z*r2inv*dU2 - (SS_d.y*SS_d.z*rinv)*r2inv*dU;
            W.fz.x = SS_d.z*SS_d.x*r2inv*dU2 - (SS_d.z*SS_d.x*rinv)*r2inv*dU;
            W.fz.y = SS_d.z*SS_d.y*r2inv*dU2 - (SS_d.z*SS_d.y*rinv)*r2inv*dU;
            W.fz.z = SS_d.z*SS_d.z*r2inv*dU2 + (SS_r - SS_d.z*SS_d.z*rinv)*r2inv*dU;



        }
    }


    return W;
}
