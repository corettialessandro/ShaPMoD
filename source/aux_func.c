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
    double matrixTimesVector[ndimension];
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

        errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

        DIRECTION_OLD[i] = RESIDUE_OLD[i];

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

            matrixTimesVector[i] = 0.;

            for (j=0; j<ndimension; j++) {

                matrixTimesVector[i] += matrix[i][j]*DIRECTION_OLD[j];

            }
        }

        for (i=0; i<ndimension; i++) {

            alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];

            alpha_denom += DIRECTION_OLD[i]*matrixTimesVector[i];

            
        }

        alpha = alpha_num/alpha_denom;


        for (i=0; i<ndimension; i++) {

            vector_x[i] += alpha*DIRECTION_OLD[i];

            RESIDUE[i] = RESIDUE_OLD[i];

            RESIDUE[i] -= (alpha*matrixTimesVector[i]);

            beta_num += RESIDUE[i]*RESIDUE[i];
            beta_denom += RESIDUE_OLD[i]*RESIDUE_OLD[i];
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
            errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

        }

        errorNorm = sqrt(errorNorm);

    }
    // printf("%.4e", errorNorm);
    
    // printf(" \n vector x = %.4e %.4e \n", vector_x[0], vector_x[1]);
    // for (i=0; i<ndimension; i++) {
    //     printf("lambda %d = %.8e\n",i, vector_x[i]);
    
    // }
    return;

}

void TrickyLinearConjugateGradient(double **Bmatrix, double *vector_b, double *vector_x, int ndimension){

    int i, j, k;
    int counterLCG = 0;
    double alpha=0., alpha_num, alpha_denom, beta=0., beta_num, beta_denom, errorNorm = 0.;
    double gammaPrime[ndimension], directionPrime[ndimension];
    double matrixTimesVector[ndimension];
    //double testb[2];


    for (i=0; i<ndimension; i++) {
        RESIDUE_OLD[i] = 0.;
        DIRECTION_OLD[i] = 0.;
        RESIDUE[i] = 0.;
        DIRECTION[i] = 0.;
        ERRORVECTOR[i] = 0.;
    }

    for (i=0; i<ndimension; i++) {

        gammaPrime[i] = 0.;

        for (j=0; j<ndimension; j++) {

            gammaPrime[i] += Bmatrix[i][j]*vector_x[j];

        }
    }


    for (i=0; i<ndimension; i++) {

        RESIDUE_OLD[i] = vector_b[i];

        for (j=0; j<ndimension; j++) {

            RESIDUE_OLD[i] -= (Bmatrix[i][j]*gammaPrime[j]);

        }

        DIRECTION_OLD[i] = RESIDUE_OLD[i];

        errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

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

            directionPrime[i] = 0.;

            for (j=0; j<ndimension; j++) {
                
                directionPrime[i] += Bmatrix[i][j]*DIRECTION_OLD[j];

            }
        }

        for (i=0; i<ndimension; i++) {

            matrixTimesVector[i] = 0.;

            for (j=0; j<ndimension; j++) {

                matrixTimesVector[i] += Bmatrix[i][j]*directionPrime[j];

            }
        }

        for (i=0; i<ndimension; i++) {

            alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];

            alpha_denom += DIRECTION_OLD[i]*matrixTimesVector[i];

        }

        alpha = alpha_num/alpha_denom;

        for (i=0; i<ndimension; i++) {

            vector_x[i] += alpha*DIRECTION_OLD[i];

            RESIDUE[i] = RESIDUE_OLD[i];

            RESIDUE[i] -= (alpha*matrixTimesVector[i]);

            
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

            errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

        }


        errorNorm = sqrt(errorNorm);



    }
    // printf("%.4e", errorNorm);
    
    // printf(" \n vector x = %.4e %.4e \n", vector_x[0], vector_x[1]);
    // for (i=0; i<ndimension; i++) {
    //     printf("lambda %d = %.8e\n",i, vector_x[i]);
    
    // }
    return;

}

void TrickyLinearConjugateGradientCellList(double **Bmatrix, double *vector_b, double *vector_x, int ndimension){

    int i, j, k;
    int counterLCG = 0;
    double alpha=0., alpha_num, alpha_denom, beta=0., beta_num, beta_denom, errorNorm = 0.;
    double gammaPrime[ndimension], directionPrime[ndimension];
    double matrixTimesVector[ndimension];
    //double testb[2];
    int neighlist[1000];
    int p;


    for (i=0; i<ndimension; i++) {
        RESIDUE_OLD[i] = 0.;
        DIRECTION_OLD[i] = 0.;
        RESIDUE[i] = 0.;
        DIRECTION[i] = 0.;
        ERRORVECTOR[i] = 0.;
        
        vector_x[i] = 0.; // seems to faster the convergence
    }

    for (i=0; i<ndimension/3.; i++) {

        gammaPrime[3*i] = 0.;
        gammaPrime[3*i+1] = 0.;
        gammaPrime[3*i+2] = 0.;

        List_Of_Neighs(i,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {
            j = neighlist[p];

            if (j < NATOMSPERSPEC[0]) {

                gammaPrime[3*i] += Bmatrix[3*i][3*j]*vector_x[3*j];
                gammaPrime[3*i] += Bmatrix[3*i][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i] += Bmatrix[3*i][3*j+2]*vector_x[3*j+2];

                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j]*vector_x[3*j];
                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j+2]*vector_x[3*j+2];

                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j]*vector_x[3*j];
                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j+2]*vector_x[3*j+2];
            }


        }
    }

    for (i=0; i<ndimension/3.; i++) {

        RESIDUE_OLD[3*i] = vector_b[3*i];
        RESIDUE_OLD[3*i+1] = vector_b[3*i+1];
        RESIDUE_OLD[3*i+2] = vector_b[3*i+2];

        List_Of_Neighs(i,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {
            j = neighlist[p];

            if (j < NATOMSPERSPEC[0]) {

                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j+2]*gammaPrime[3*j+2];

                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j+2]*gammaPrime[3*j+2];

                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j+2]*gammaPrime[3*j+2];
            }

        }

        DIRECTION_OLD[3*i] = RESIDUE_OLD[3*i];
        DIRECTION_OLD[3*i+1] = RESIDUE_OLD[3*i+1];
        DIRECTION_OLD[3*i+2] = RESIDUE_OLD[3*i+2];

        errorNorm += RESIDUE_OLD[3*i]*RESIDUE_OLD[3*i] + RESIDUE_OLD[3*i+1]*RESIDUE_OLD[3*i+1] + RESIDUE_OLD[3*i+2]*RESIDUE_OLD[3*i+2];
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


        for (i=0; i<ndimension/3.; i++) {

            directionPrime[3*i] = 0.;
            directionPrime[3*i+1] = 0.;
            directionPrime[3*i+2] = 0.;

            List_Of_Neighs(i,neighlist,1);
            for (p=1;p<=neighlist[0];p++) {
                j = neighlist[p];

                if (j < NATOMSPERSPEC[0]) {

                    directionPrime[3*i] += Bmatrix[3*i][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i] += Bmatrix[3*i][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i] += Bmatrix[3*i][3*j+2]*DIRECTION_OLD[3*j+2];

                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j+2]*DIRECTION_OLD[3*j+2];

                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j+2]*DIRECTION_OLD[3*j+2];

                }


            }
        }


        for (i=0; i<ndimension/3.; i++) {

            matrixTimesVector[3*i] = 0.;
            matrixTimesVector[3*i+1] = 0.;
            matrixTimesVector[3*i+2] = 0.;

            List_Of_Neighs(i,neighlist,1);
            for (p=1;p<=neighlist[0];p++) {
                j = neighlist[p];

                if (j < NATOMSPERSPEC[0]) {

                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j+2]*directionPrime[3*j+2];

                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j+2]*directionPrime[3*j+2];

                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j+2]*directionPrime[3*j+2];
                }

            }
        }

        for (i=0; i<ndimension; i++) {

            alpha_num += RESIDUE_OLD[i]*RESIDUE_OLD[i];

            alpha_denom += DIRECTION_OLD[i]*matrixTimesVector[i];

        }

        alpha = alpha_num/alpha_denom;

        for (i=0; i<ndimension; i++) {

            vector_x[i] += alpha*DIRECTION_OLD[i];

            RESIDUE[i] = RESIDUE_OLD[i];

            RESIDUE[i] -= (alpha*matrixTimesVector[i]);

            
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

            errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

        }


        errorNorm = sqrt(errorNorm);



    }
    printf("%.d \n", counterLCG);
    
    // printf(" \n vector x = %.4e %.4e \n", vector_x[0], vector_x[1]);
    // for (i=0; i<ndimension; i++) {
    //     printf("lambda %d = %.8e\n",i, vector_x[i]);
    
    // }
    return;

}

void TrickyPreconditionedLinearConjugateGradientCellList(double **Bmatrix, double **Cmatrix, double *vector_b, double *vector_x, int ndimension){

    int i, j, k;
    int counterLCG = 0;
    double alpha=0., alpha_num, alpha_denom, beta=0., beta_num, beta_denom, errorNorm = 0.;
    double gammaPrime[ndimension], directionPrime[ndimension], preconditioned[ndimension], preconditioned_old[ndimension];
    double matrixTimesVector[ndimension];
    //double testb[2];
    int neighlist[1000];
    int p;


    for (i=0; i<ndimension; i++) {
        RESIDUE_OLD[i] = 0.;
        DIRECTION_OLD[i] = 0.;
        RESIDUE[i] = 0.;
        DIRECTION[i] = 0.;
        ERRORVECTOR[i] = 0.;
        preconditioned[i] = 0.;
        preconditioned_old[i]= 0.;
        vector_x[i] = 0.; // seems to faster the convergence
    }

    for (i=0; i<ndimension/3.; i++) {

        gammaPrime[3*i] = 0.;
        gammaPrime[3*i+1] = 0.;
        gammaPrime[3*i+2] = 0.;

        List_Of_Neighs(i,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {
            j = neighlist[p];

            if (j < NATOMSPERSPEC[0]) {

                gammaPrime[3*i] += Bmatrix[3*i][3*j]*vector_x[3*j];
                gammaPrime[3*i] += Bmatrix[3*i][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i] += Bmatrix[3*i][3*j+2]*vector_x[3*j+2];

                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j]*vector_x[3*j];
                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i+1] += Bmatrix[3*i+1][3*j+2]*vector_x[3*j+2];

                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j]*vector_x[3*j];
                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j+1]*vector_x[3*j+1];
                gammaPrime[3*i+2] += Bmatrix[3*i+2][3*j+2]*vector_x[3*j+2];
            }


        }
    }

    for (i=0; i<ndimension/3.; i++) {

        RESIDUE_OLD[3*i] = vector_b[3*i];
        RESIDUE_OLD[3*i+1] = vector_b[3*i+1];
        RESIDUE_OLD[3*i+2] = vector_b[3*i+2];

        List_Of_Neighs(i,neighlist,1);
        for (p=1;p<=neighlist[0];p++) {
            j = neighlist[p];

            if (j < NATOMSPERSPEC[0]) {

                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i] -= Bmatrix[3*i][3*j+2]*gammaPrime[3*j+2];

                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i+1] -= Bmatrix[3*i+1][3*j+2]*gammaPrime[3*j+2];

                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j]*gammaPrime[3*j];
                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j+1]*gammaPrime[3*j+1];
                RESIDUE_OLD[3*i+2] -= Bmatrix[3*i+2][3*j+2]*gammaPrime[3*j+2];
            }

        }

        preconditioned_old[3*i] = 1/Cmatrix[3*i][3*i]*RESIDUE_OLD[3*i];
        preconditioned_old[3*i + 1] = 1/Cmatrix[3*i + 1][3*i + 1]*RESIDUE_OLD[3*i + 1];
        preconditioned_old[3*i + 2] = 1/Cmatrix[3*i + 2][3*i + 2]*RESIDUE_OLD[3*i + 2];

        DIRECTION_OLD[3*i] = RESIDUE_OLD[3*i];
        DIRECTION_OLD[3*i+1] = RESIDUE_OLD[3*i+1];
        DIRECTION_OLD[3*i+2] = RESIDUE_OLD[3*i+2];
        

        errorNorm += RESIDUE_OLD[3*i]*RESIDUE_OLD[3*i] + RESIDUE_OLD[3*i+1]*RESIDUE_OLD[3*i+1] + RESIDUE_OLD[3*i+2]*RESIDUE_OLD[3*i+2];
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


        for (i=0; i<ndimension/3.; i++) {

            directionPrime[3*i] = 0.;
            directionPrime[3*i+1] = 0.;
            directionPrime[3*i+2] = 0.;

            List_Of_Neighs(i,neighlist,1);
            for (p=1;p<=neighlist[0];p++) {
                j = neighlist[p];

                if (j < NATOMSPERSPEC[0]) {

                    directionPrime[3*i] += Bmatrix[3*i][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i] += Bmatrix[3*i][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i] += Bmatrix[3*i][3*j+2]*DIRECTION_OLD[3*j+2];

                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i+1] += Bmatrix[3*i+1][3*j+2]*DIRECTION_OLD[3*j+2];

                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j]*DIRECTION_OLD[3*j];
                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j+1]*DIRECTION_OLD[3*j+1];
                    directionPrime[3*i+2] += Bmatrix[3*i+2][3*j+2]*DIRECTION_OLD[3*j+2];

                }


            }
        }


        for (i=0; i<ndimension/3.; i++) {

            matrixTimesVector[3*i] = 0.;
            matrixTimesVector[3*i+1] = 0.;
            matrixTimesVector[3*i+2] = 0.;

            List_Of_Neighs(i,neighlist,1);
            for (p=1;p<=neighlist[0];p++) {
                j = neighlist[p];

                if (j < NATOMSPERSPEC[0]) {

                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i] += Bmatrix[3*i][3*j+2]*directionPrime[3*j+2];

                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i+1] += Bmatrix[3*i+1][3*j+2]*directionPrime[3*j+2];

                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j]*directionPrime[3*j];
                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j+1]*directionPrime[3*j+1];
                    matrixTimesVector[3*i+2] += Bmatrix[3*i+2][3*j+2]*directionPrime[3*j+2];
                }

            }
        }

        for (i=0; i<ndimension; i++) {

            alpha_num += RESIDUE_OLD[i]*preconditioned_old[i];

            alpha_denom += DIRECTION_OLD[i]*matrixTimesVector[i];

        }

        alpha = alpha_num/alpha_denom;

        for (i=0; i<ndimension; i++) {

            vector_x[i] += alpha*DIRECTION_OLD[i];

            RESIDUE[i] = RESIDUE_OLD[i];

            RESIDUE[i] -= (alpha*matrixTimesVector[i]);

            preconditioned[i] = 1/Cmatrix[i][i]*RESIDUE[i];
            //printf("invC = %.4e \n", 1/Cmatrix[i][i]);

        }

        for (i=0; i<ndimension; i++) {

            beta_num += RESIDUE[i]*preconditioned[i];
            beta_denom += RESIDUE_OLD[i]*preconditioned_old[i];


        }
        beta = beta_num/beta_denom;
        


        for (i=0; i<ndimension; i++) {

            DIRECTION[i] = preconditioned[i] + beta*DIRECTION_OLD[i];

        }

        for (i=0; i<ndimension; i++) {

            DIRECTION_OLD[i] = DIRECTION[i];
            RESIDUE_OLD[i] = RESIDUE[i];

            errorNorm += RESIDUE_OLD[i]*RESIDUE_OLD[i];

        }


        errorNorm = sqrt(errorNorm);
        printf("error = %.4e \n", errorNorm);



    }
    printf("%.d \n", counterLCG);
    exit(0);
    // printf(" \n vector x = %.4e %.4e \n", vector_x[0], vector_x[1]);
    // for (i=0; i<ndimension; i++) {
    //     printf("lambda %d = %.8e\n",i, vector_x[i]);
    
    // }
    return;

}

