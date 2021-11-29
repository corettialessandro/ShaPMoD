//
//  aux_func.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "aux_func.h"

struct point Distance(struct point A, struct point B){

    struct point AB_d;

    double AB_r2;

    AB_d.x = A.x - B.x;
    AB_d.y = A.y - B.y;
    AB_d.z = A.z - B.z;

    if (R_CUT > 0) {

        AB_d.x -= LBOX*nearbyint(AB_d.x/LBOX);
        AB_d.y -= LBOX*nearbyint(AB_d.y/LBOX);
        AB_d.z -= LBOX*nearbyint(AB_d.z/LBOX);

        if ((AB_r2 = modsq(AB_d)) > 3./4.*LBOX*LBOX) {

//            printf("\naux_func.c -> Distance() ERROR: Distance greater than half box length.\nLBOX/2 = %lf\td = (%lf, %lf, %lf)\nExecution aborted.\n\n", .5*LBOX, AB_d.x, AB_d.y, AB_d.z);
//            exit(EXIT_FAILURE);
        }
    }

    return AB_d;
}

struct point d_rirj(struct point ri, struct point rj){

    struct point d;

    d = Distance(ri, rj);

    return d;
}

struct point d_rirhoj(struct point ri, struct point rhoj, struct point rj){

    struct point d, sj;

    d = Distance(ri, rj);
    sj = Distance(rhoj, rj);

    d.x -= sj.x;
    d.y -= sj.y;
    d.z -= sj.z;

    return d;
}

struct point d_rhoirj(struct point rhoi, struct point ri, struct point rj){

    struct point d, si;

    d = Distance(ri, rj);
    si = Distance(rhoi, ri);

    d.x += si.x;
    d.y += si.y;
    d.z += si.z;

    return d;
}

struct point d_rhoirhoj(struct point rhoi, struct point ri, struct point rhoj, struct point rj){

    struct point d, si, sj;

    d = Distance(ri, rj);
    si = Distance(rhoi, ri);
    sj = Distance(rhoj, rj);

    d.x += si.x - sj.x;
    d.y += si.y - sj.y;
    d.z += si.z - sj.z;

    return d;
}

double Inner(struct point A, struct point B){

    return (A.x*B.x + A.y*B.y + A.z*B.z);
}

struct point Vector(struct point A, struct point B){

    struct point C;

    C.x = A.y*B.z - A.z*B.y;
    C.y = A.z*B.x - A.x*B.z;
    C.z = A.x*B.y - A.y*B.x;

    return C;
}

double mod(struct point A){

    return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

double modsq(struct point A){

    return (A.x*A.x + A.y*A.y + A.z*A.z);
}

struct point RNDM_Rotate(struct point r){

    struct point r_rot;

    double theta, phi;

    theta = lrand48()/(RAND_MAX+1.)*M_PI;
    phi = lrand48()/(RAND_MAX+1.)*2.*M_PI;

    r_rot.x = r.x*sin(theta)*cos(phi);
    r_rot.y = r.y*sin(theta)*sin(phi);
    r_rot.z = r.z*cos(theta);

    return r_rot;
}

double Variance(struct point v[], int n){

    int i;
    double var = 0.0;

    for (i=0; i<n; i++) {

        var += (v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z);
    }

    return sqrt(var/n);
}

double OTRAvg(int N, double Estimator, double Avg){

    return Avg = (Avg*(N-1) + Estimator)/(double)N;
}

double Thermostat(struct point r_tm1[], struct point r_t[],  double temp){

    int i;
    double v2 = 0, alpha, g = 3.*(NPART - 1.);
    struct point vi;

    for (i=0; i<NPART; i++) {

        vi.x = (r_t[i].x - r_tm1[i].x)/DT;
        vi.y = (r_t[i].y - r_tm1[i].y)/DT;
        vi.z = (r_t[i].z - r_tm1[i].z)/DT;

        v2 += M[INDX[i]]*modsq(vi);
    }

    return alpha = sqrt((1.-_TEMP_TOL)*temp*g/v2);
}

void LinearConjugateGradient(double **matrix, double *vector_b, double *vector_x, int ndimension){

    int i, j, k;
    int counterLCG = 0;
    double alpha=0., alpha_num, alpha_denom, beta=0., beta_num, beta_denom, errorNorm = 0.;
    //double testb[2];


    for (i=0; i<ndimension; i++) {
        RESIDUE_OLD[i] = 0.;
        DIRECTION_OLD[i] = 0.;
        RESIDUE[i] = 0.;
        DIRECTION[i] = 0.;
        ERRORVECTOR[i] = 0.;
    }

    for (i=0; i<ndimension; i++) {

        RESIDUE_OLD[i] = vector_b[i];

        for (j=0; j<ndimension; j++) {

            RESIDUE_OLD[i] -= (matrix[i][j]*vector_x[j]);

        }

        DIRECTION_OLD[i] = RESIDUE_OLD[i];

    }


    for (i=0; i<ndimension; i++) {

        ERRORVECTOR[i] = vector_b[i];

        for (j=0; j<ndimension; j++) {

            ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j]);

        }
        errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
    }

    errorNorm = sqrt(errorNorm);

    while (errorNorm > LCG_TOL) {

        alpha_num = 0.;
        alpha_denom = 0.;
        beta_num = 0.;
        beta_denom = 0.;
        errorNorm = 0.;

        counterLCG++;

        if (counterLCG > _MAX_ITER) {

            printf("\n Iteration limit exceeded for LinearConjugateGradient \n");
            exit(EXIT_FAILURE);

        }

        for (i=0; i<ndimension; i++) {

            alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];

            for (j=0; j<ndimension; j++) {

                alpha_denom += DIRECTION_OLD[i]*matrix[i][j]*DIRECTION_OLD[j];

            }
        }

        alpha = alpha_num/alpha_denom;


        for (i=0; i<ndimension; i++) {

            vector_x[i] += alpha*DIRECTION_OLD[i];
            RESIDUE[i] = RESIDUE_OLD[i];

            for (j=0; j<ndimension; j++) {

                RESIDUE[i] -= (alpha*matrix[i][j]*DIRECTION_OLD[j]);

            }
        }

        for (i=0; i<ndimension; i++) {

            beta_num += RESIDUE[i]*RESIDUE[i];
            beta_denom += RESIDUE_OLD[i]*RESIDUE_OLD[i];


        }
        beta = beta_num/beta_denom;

        for (i=0; i<ndimension; i++) {

            DIRECTION[i] = RESIDUE[i] + beta*DIRECTION_OLD[i];

        }

        for (i=0; i<ndimension; i++) {

            DIRECTION_OLD[i] = DIRECTION[i];
            RESIDUE_OLD[i] = RESIDUE[i];

        }

        for (i=0; i<ndimension; i++) {

            ERRORVECTOR[i] = vector_b[i];

            for (j=0; j<ndimension; j++) {

                ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j]);

            }
            errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
        }

        errorNorm = sqrt(errorNorm);



    }
    // printf("%.4e", errorNorm);
    //
    //printf(" \n vector x = %.4e %.4e \n", vector_x[0], vector_x[1]);
    // for (i=0; i<ndimension; i++) {
    //     printf("lambda %d = %.8e\n",i, vector_x[i]);
    //
    // }
    return;
// ------------------- y ----------------------
//     if (component == 'y'){
//
//         for (i=0; i<ndimension; i++) {
//             RESIDUE_OLD[i] = 0.;
//             DIRECTION_OLD[i] = 0.;
//             RESIDUE[i] = 0.;
//             DIRECTION[i] = 0.;
//             ERRORVECTOR[i] = 0.;
//         }
//
//         for (i=0; i<ndimension; i++) {
//
//             RESIDUE_OLD[i] = vector_b[i].y;
//
//             for (j=0; j<ndimension; j++) {
//
//                 RESIDUE_OLD[i] -= (matrix[i][j]*vector_x[j].y);
//
//             }
//
//             DIRECTION_OLD[i] = RESIDUE_OLD[i];
//
//         }
//
//
//         for (i=0; i<ndimension; i++) {
//
//             ERRORVECTOR[i] = vector_b[i].y;
//
//             for (j=0; j<ndimension; j++) {
//
//                 ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j].y);
//
//             }
//             errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
//         }
//
//         errorNorm = sqrt(errorNorm);
//
//         while (errorNorm > LCG_TOL) {
//
//             alpha_num = 0.;
//             alpha_denom = 0.;
//             beta_num = 0.;
//             beta_denom = 0.;
//             errorNorm = 0.;
//
//             counterLCG++;
//
//             if (counterLCG > _MAX_ITER) {
//
//                 printf("\n Iteration limit exceeded for LinearConjugateGradient \n");
//                 exit(EXIT_FAILURE);
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];
//
//                 for (j=0; j<ndimension; j++) {
//
//                     alpha_denom += DIRECTION_OLD[i]*matrix[i][j]*DIRECTION_OLD[j];
//
//                 }
//             }
//
//             alpha = alpha_num/alpha_denom;
//
//
//             for (i=0; i<ndimension; i++) {
//
//                 vector_x[i].y += alpha*DIRECTION_OLD[i];
//                 RESIDUE[i] = RESIDUE_OLD[i];
//
//                 for (j=0; j<ndimension; j++) {
//
//                     RESIDUE[i] -= (alpha*matrix[i][j]*DIRECTION_OLD[j]);
//
//                 }
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 beta_num += RESIDUE[i]*RESIDUE[i];
//                 beta_denom += RESIDUE_OLD[i]*RESIDUE_OLD[i];
//
//
//             }
//             beta = beta_num/beta_denom;
//
//             for (i=0; i<ndimension; i++) {
//
//                 DIRECTION[i] = RESIDUE[i] + beta*DIRECTION_OLD[i];
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 DIRECTION_OLD[i] = DIRECTION[i];
//                 RESIDUE_OLD[i] = RESIDUE[i];
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 ERRORVECTOR[i] = vector_b[i].y;
//
//                 for (j=0; j<ndimension; j++) {
//
//                     ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j].y);
//
//                 }
//                 errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
//             }
//
//             errorNorm = sqrt(errorNorm);
//
//
//
//         }
//         //printf(" \n vector x = %.4e %.4e \n", vector_x[0].y, vector_x[1].y);
//         return;
//     }
// // ------------------- z ----------------------
//     if (component == 'z'){
//
//         for (i=0; i<ndimension; i++) {
//             RESIDUE_OLD[i] = 0.;
//             DIRECTION_OLD[i] = 0.;
//             RESIDUE[i] = 0.;
//             DIRECTION[i] = 0.;
//             ERRORVECTOR[i] = 0.;
//         }
//
//         for (i=0; i<ndimension; i++) {
//
//             RESIDUE_OLD[i] = vector_b[i].z;
//
//             for (j=0; j<ndimension; j++) {
//
//                 RESIDUE_OLD[i] -= (matrix[i][j]*vector_x[j].z);
//
//             }
//
//             DIRECTION_OLD[i] = RESIDUE_OLD[i];
//
//         }
//
//
//         for (i=0; i<ndimension; i++) {
//
//             ERRORVECTOR[i] = vector_b[i].z;
//
//             for (j=0; j<ndimension; j++) {
//
//                 ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j].z);
//
//             }
//             errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
//         }
//
//         errorNorm = sqrt(errorNorm);
//
//         while (errorNorm > LCG_TOL) {
//
//             alpha_num = 0.;
//             alpha_denom = 0.;
//             beta_num = 0.;
//             beta_denom = 0.;
//             errorNorm = 0.;
//
//             counterLCG++;
//
//             if (counterLCG > _MAX_ITER) {
//
//                 printf("\n Iteration limit exceeded for LinearConjugateGradient \n");
//                 exit(EXIT_FAILURE);
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];
//
//                 for (j=0; j<ndimension; j++) {
//
//                     alpha_denom += DIRECTION_OLD[i]*matrix[i][j]*DIRECTION_OLD[j];
//
//                 }
//             }
//
//             alpha = alpha_num/alpha_denom;
//
//
//             for (i=0; i<ndimension; i++) {
//
//                 vector_x[i].z += alpha*DIRECTION_OLD[i];
//                 RESIDUE[i] = RESIDUE_OLD[i];
//
//                 for (j=0; j<ndimension; j++) {
//
//                     RESIDUE[i] -= (alpha*matrix[i][j]*DIRECTION_OLD[j]);
//
//                 }
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 beta_num += RESIDUE[i]*RESIDUE[i];
//                 beta_denom += RESIDUE_OLD[i]*RESIDUE_OLD[i];
//
//
//             }
//             beta = beta_num/beta_denom;
//
//             for (i=0; i<ndimension; i++) {
//
//                 DIRECTION[i] = RESIDUE[i] + beta*DIRECTION_OLD[i];
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 DIRECTION_OLD[i] = DIRECTION[i];
//                 RESIDUE_OLD[i] = RESIDUE[i];
//
//             }
//
//             for (i=0; i<ndimension; i++) {
//
//                 ERRORVECTOR[i] = vector_b[i].z;
//
//                 for (j=0; j<ndimension; j++) {
//
//                     ERRORVECTOR[i] -= (matrix[i][j]*vector_x[j].z);
//
//                 }
//                 errorNorm += ERRORVECTOR[i]*ERRORVECTOR[i];
//             }
//
//             errorNorm = sqrt(errorNorm);
//
//
//
//         }
//         //printf(" \n vector x = %.4e %.4e \n", vector_x[0].z, vector_x[1].z);
//         return;
//     }
//


}
