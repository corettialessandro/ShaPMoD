//
//  common.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "common.h"

struct point *PARTPOS_TM1;
struct point *PARTPOS_T;
struct point *PARTPOS_TP1;
struct point *PARTVEL;
struct point *PARTMOM_T;
struct point *PARTMOM_TP05;
struct point *PARTMOM_TP1;
struct point *SHELLPOS_TM1;
struct point *SHELLPOS_T;
struct point *SHELLPOS_TP1;
struct point *SHELLVEL;
struct point *SHELLVEL_TP1;
struct point *SHELLACC_TM1;

struct tensor **DPHIDRHO_T;
struct tensor **DPHIDVRHO_T;
struct point *GAMMA;

struct point *PHI;
struct point *PHI_OLD;
struct point *SEARCHDIR;
struct point *RHO_OLD;

char MODE;
char SRMODE;
char EWALD;
char POT;

int NMOLECULES;
int NATOMSPEC;
int *NATOMSPERSPEC;
int NATOMSPERMOLECULE;
int NPART;
int NINTER;

double DT;
int NTIMESTEPS;

double DENSITY;
double *DENS_RED;
double TEMP;

double R_CUT;
double LBOX;
double RCUT;
double RCUT_SQ;
double LR_RCUT;
double ALPHA;

double B;
double LAMBDA;

int *INDX;

char **NAME;
double *M;
double *MU;
double *Q;
double *CHI;
double *K;
double *SIGMA;

double MTOT = 0;

double *A;
double *C;
double *D;

double LR_VCUT;
double EW_VCUT;
double *C_VCUT;
double *S_VCUT;
double LR_DVCUT;
double EW_DVCUT;
double *C_DVCUT;
double *S_DVCUT;
double LR_DDVCUT;
double EW_DDVCUT;
double *C_DDVCUT;
double *S_DDVCUT;
double *C_VTAIL;
double *S_VTAIL;
double *C_PTAIL;
double *S_PTAIL;

int **GCOUNT;
int NGCOUNT;
double **GOFR;
double DG;

double LOW_TOL;
int SR_ITERS = 0;
double SR_DISCR = 0;
int GET_OUT = 0;

char RUNNAME[_MAX_STR_LENGTH] = "";
char INPUTFOL[_MAX_STR_LENGTH] = "input/default/";
char INPUTFN[_MAX_STR_LENGTH] = "input/default/input.txt";
char OUTPUTFOL[_MAX_STR_LENGTH] = "output/default/";
char PRECONFFN[_MAX_STR_LENGTH] = "input/default/checkpoint.txt";
char INITCONFFN[_MAX_STR_LENGTH] = "input/default/init_checkpoint.txt";
char SAVECONFFN[_MAX_STR_LENGTH] = "input/default/save_conf.txt";
char SAVECONFFNM1[_MAX_STR_LENGTH] = "input/default/save_confm1.txt";

int IANSCREEN;
int IANFILE;
int IPS;
int IVMD;
int ISHAKE;
int IGOFR;
int ICHECK;

int DEBUG_FLAG = 0;
int VERBOSE_FLAG = 0;
int OUTPUT_FLAG = 0;
int ANALYSIS_FLAG = 0;
int PRECONFIG_FLAG = 0;
int COLLIDING_FLAG = 0;

char THERMOSTAT;

double ITEMP = 0;

int NANCOUNT = 0;
double ETOTAVG = 0;
double EKINAVG = 0;
double EPOTAVG = 0;
double EPOLAVG = 0;
double TEMPAVG = 0;
double PRESSAVG = 0;
double ETOTAVGSQ = 0;
double EKINAVGSQ = 0;
double EPOTAVGSQ = 0;
double EPOLAVGSQ = 0;
double TEMPAVGSQ = 0;
double PRESSAVGSQ = 0;
struct point ANMOMAVG = {0};

struct point *CF;
struct point *SF;

double EPOT = 0;
double PRESS = 0;

void ReadInput(){

    FILE *fp_input;

    int i, j, k;
    int n1, n2;
    int pointer_flag = 0, tnpart, *dens_red_count, ilr_vcut = 1., isr_vcut = 1., i_vtail = 1.;
    double tlbox, sigmaij;
    char dummy;

    if ((fp_input = fopen(INPUTFN, "r+")) == NULL){
        printf("\ninput.c -> ReadInput() Error: File '%s' not found!\n", INPUTFN);
        exit(EXIT_FAILURE);
    }

    fscanf(fp_input, "%c %*[^\n]\n", &MODE);
    fscanf(fp_input, "%c %*[^\n]\n", &SRMODE);
    fscanf(fp_input, "%c %*[^\n]\n", &EWALD);
    fscanf(fp_input, "%c %*[^\n]\n", &POT);
    fscanf(fp_input, "%d %*[^\n]", &NMOLECULES);
    fscanf(fp_input, "%d %*[^\n]", &NATOMSPEC);

    if ((NATOMSPERSPEC = (int *) calloc(NATOMSPEC, sizeof(int))) == NULL) pointer_flag = 1;
    NATOMSPERMOLECULE = 0;

    for (i=0; i<NATOMSPEC; i++) {

        fscanf(fp_input, "%d %*[^\n]", &NATOMSPERSPEC[i]);
        NATOMSPERMOLECULE += NATOMSPERSPEC[i];
    }

    NPART = NATOMSPERMOLECULE*NMOLECULES;
    NINTER = (int)round(NATOMSPEC*(NATOMSPEC+1.)/2.);
    // allocate memory for arrays and check if done properly
    if ((PARTPOS_TM1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 2;
    if ((PARTPOS_T = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 3;
    if ((PARTPOS_TP1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 4;
    if ((PARTMOM_T = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 5;
    if ((PARTMOM_TP05 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 5;
    if ((PARTMOM_TP1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 5;
    if ((PARTVEL = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 5;
    if ((SHELLPOS_TM1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 6;
    if ((SHELLPOS_T = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 7;
    if ((SHELLPOS_TP1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 8;
    if ((SHELLVEL = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 8;
    if ((SHELLVEL_TP1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 8;
    if ((SHELLACC_TM1 = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 8;

    if (EWALD == 'T') {
        if ((CF = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 80;
        if ((SF = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 81;
    }

    if (SRMODE == 'S') {

        if ((DPHIDRHO_T = (struct tensor **)calloc(NPART, sizeof(struct tensor *))) == NULL) pointer_flag = 100;
        if ((DPHIDVRHO_T = (struct tensor **)calloc(NPART, sizeof(struct tensor *))) == NULL) pointer_flag = 100;
        if ((GAMMA = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 9;

        for (i=0; i<NPART; i++) {

            if ((DPHIDRHO_T[i] = (struct tensor *)calloc(NPART, sizeof(struct tensor))) == NULL) pointer_flag = 100+1+i;
            if ((DPHIDVRHO_T[i] = (struct tensor *)calloc(NPART, sizeof(struct tensor))) == NULL) pointer_flag = 100+1+i;
        }
    }

    if (SRMODE == 'C') {

        if ((PHI = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 2;
        if ((PHI_OLD = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 2;
        if ((SEARCHDIR = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 2;
        if ((RHO_OLD = (struct point *)calloc(NPART, sizeof(struct point))) == NULL) pointer_flag = 2;
    }

    fscanf(fp_input, "%lf %*[^\n]", &DT);
    DT /= _TIME_CONV;
    fscanf(fp_input, "%d %*[^\n]", &NTIMESTEPS);

    fscanf(fp_input, "%lf %*[^\n]\n", &DENSITY);
    DENSITY /= _DENS_CONV;
    fscanf(fp_input, "%c %*[^\n]", &THERMOSTAT);
    fscanf(fp_input, "%lf %*[^\n]", &TEMP);
    TEMP /= _TEMP_CONV;
    fscanf(fp_input, "%lf %*[^\n]", &R_CUT);

    LBOX = pow(((double)NPART/DENSITY), 1./3.);
    RCUT = R_CUT*LBOX;
    LR_RCUT = .5*LBOX;

    //R_CUT < 0 && R_CUT != -1 -> NO RCUT
    if (R_CUT <= 0) {

        RCUT = .5*LBOX;

        //R_CUT == -1 -> NO PBC, NO RCUT
        if (R_CUT == -1) ilr_vcut = isr_vcut = i_vtail = 0;
    }

    RCUT_SQ = RCUT*RCUT;

    if (EWALD == 'T') ALPHA = _ALPHA/LBOX;

    if ((NAME = (char **)calloc(NATOMSPEC, sizeof(char *))) == NULL) pointer_flag = 1000;
    if ((M = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 10;
    if ((MU = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 11;
    if ((Q = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 12;
    if ((CHI = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 13;
    if ((K = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 14;
    if ((SIGMA = (double *)calloc(NATOMSPEC, sizeof(double))) == NULL) pointer_flag = 15;

    for (i=0; i<NATOMSPEC; i++) {

        if ((NAME[i] = (char *)calloc(3, sizeof(char))) == NULL) pointer_flag = 1000+1+i;

        fscanf(fp_input, "%s %*[^\n]", NAME[i]);
        fscanf(fp_input, "%lf %*[^\n]", &M[i]);
        fscanf(fp_input, "%lf %*[^\n]", &MU[i]);
        fscanf(fp_input, "%lf %*[^\n]", &Q[i]);
        fscanf(fp_input, "%lf %*[^\n]", &CHI[i]);

        if (MODE == 'S') {

            Q[i] += CHI[i];
            CHI[i] = MU[i] = K[i] = 0.;
        }

        fscanf(fp_input, "%lf %*[^\n]", &SIGMA[i]);
        fscanf(fp_input, "%lf %*[^\n]", &K[i]);
    }

    if ((INDX = (int *)calloc(NPART, sizeof(int))) == NULL) pointer_flag = 16;

    for (i=0; i<NPART-(NATOMSPERMOLECULE-1); i+=NATOMSPERMOLECULE) {

        for (j=0; j<NATOMSPEC; j++) {

            for (k=0; k<NATOMSPERSPEC[j]; k++) {

                INDX[i+j+k] = j;
                MTOT += M[INDX[i+j+k]];
            }
        }
    }

    fscanf(fp_input, "%lf %*[^\n]", &B);
    fscanf(fp_input, "%lf %*[^\n]", &LAMBDA);

    if ((A = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 17;
    if ((C = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 18;
    if ((D = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 19;

    if ((GCOUNT = (int **)calloc(NINTER, sizeof(int *))) == NULL) pointer_flag = 2000;
    if ((GOFR = (double **)calloc(NINTER, sizeof(double *))) == NULL) pointer_flag = 2100;
    DG = LBOX/(2.*_NHIST);

    if ((C_VCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((S_VCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((C_DVCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((S_DVCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((C_DDVCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((S_DDVCUT = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 21;
    if ((C_VTAIL = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 25;
    if ((S_VTAIL = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 26;
    if ((C_PTAIL = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 27;
    if ((S_PTAIL = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 28;

    for (i=0; i<NINTER; i++) {

        fscanf(fp_input, "%lf %*[^\n]", &A[i]);
        fscanf(fp_input, "%lf %*[^\n]", &C[i]);
        fscanf(fp_input, "%lf %*[^\n]", &D[i]);

        if ((GCOUNT[i] = (int *)calloc(_NHIST, sizeof(int))) == NULL) pointer_flag = 2000+1+i;
        if ((GOFR[i] = (double *)calloc(_NHIST, sizeof(double))) == NULL) pointer_flag = 2100+1+i;
    }

    if ((DENS_RED = (double *)calloc(NINTER, sizeof(double))) == NULL) pointer_flag = 29;
    if ((dens_red_count = (int *)calloc(NINTER, sizeof(int))) == NULL) pointer_flag = 30;

    for (i=0; i<NATOMSPEC; i++) {
        for (j=0; j<NATOMSPEC; j++) {

            dens_red_count[i+j]++;
            DENS_RED[i+j] = dens_red_count[i+j]*(double)NATOMSPERSPEC[i]/(double)NATOMSPERMOLECULE*(double)NATOMSPERSPEC[j]/(double)NATOMSPERMOLECULE;
        }
    }

    free(dens_red_count);

    for (i=0; i<NINTER; i++) {

        n1 = (int) floor(i/2.);
        n2 = (int) floor((i+1)/2.);

        sigmaij = (SIGMA[n1] + SIGMA[n2]);

        A[i] = B*A[i]*exp(sigmaij/LAMBDA);

        LR_VCUT = ilr_vcut/LR_RCUT;
        LR_DVCUT = -LR_VCUT/LR_RCUT;
        LR_DDVCUT = -2*LR_DVCUT/LR_RCUT;

        EW_VCUT = erfc(ALPHA*RCUT)/RCUT;
        EW_DVCUT = -M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*RCUT_SQ)/RCUT - erfc(ALPHA*RCUT)/RCUT_SQ;
        EW_DDVCUT = 2*M_2_SQRTPI*ALPHA*exp(-ALPHA*ALPHA*RCUT_SQ)/RCUT_SQ + 2*M_2_SQRTPI*ALPHA*ALPHA*ALPHA*exp(-ALPHA*ALPHA*RCUT_SQ) + 2*erfc(ALPHA*RCUT)/RCUT_SQ/RCUT;

        if (POT == 'J') {

            C_VCUT[i] = isr_vcut*(-C[i]/pow(RCUT, 6) - D[i]/pow(RCUT, 8));
            C_DVCUT[i] = isr_vcut*(6*C[i]/pow(RCUT, 7) + 8*D[i]/pow(RCUT, 9));
            C_DDVCUT[i] = isr_vcut*(-42*C[i]/pow(RCUT, 8) - 72*D[i]/pow(RCUT, 10));

            S_VCUT[i] = isr_vcut*(A[i]*exp(-RCUT/LAMBDA));
            S_DVCUT[i] = isr_vcut*(-A[i]/LAMBDA*exp(-RCUT/LAMBDA));
            S_DDVCUT[i] = isr_vcut*(A[i]/LAMBDA/LAMBDA*exp(-RCUT/LAMBDA));

            C_VTAIL[i] = i_vtail*(2*M_PI*DENS_RED[i]*NPART*DENSITY*(-C[i]/3./pow(RCUT, 3) - D[i]/5./pow(RCUT, 5)));
            S_VTAIL[i] = i_vtail*(2*M_PI*DENS_RED[i]*NPART*DENSITY*(A[i]*LAMBDA*(RCUT_SQ + 2*RCUT*LAMBDA + 2*LAMBDA*LAMBDA)*exp(-RCUT/LAMBDA)));

            C_PTAIL[i] = i_vtail*(-4./3.*M_PI*DENS_RED[i]*DENSITY*DENSITY/RCUT/RCUT/RCUT*(C[i] + 4.*D[i]/5./RCUT/RCUT));
            S_PTAIL[i] = i_vtail*(2./3.*M_PI*DENS_RED[i]*DENSITY*DENSITY*A[i]*exp(-RCUT/LAMBDA)*(RCUT*RCUT*RCUT + 3.*RCUT*RCUT*LAMBDA + 6.*RCUT*LAMBDA*LAMBDA + 6.*LAMBDA*LAMBDA*LAMBDA));

        } else if (POT == 'C'){

            C_VCUT[i] = isr_vcut*(A[i]*exp(-RCUT/LAMBDA) - C[i]/pow(RCUT, 6) - D[i]/pow(RCUT, 8));
            C_DVCUT[i] = isr_vcut*(-A[i]/LAMBDA*exp(-RCUT/LAMBDA) + 6*C[i]/pow(RCUT, 7) + 8*D[i]/pow(RCUT, 9));
            C_DDVCUT[i] = isr_vcut*(A[i]/LAMBDA/LAMBDA*exp(-RCUT/LAMBDA) - 42*C[i]/pow(RCUT, 8) - 72*D[i]/pow(RCUT, 10));

            S_VCUT[i] = 0;
            S_DVCUT[i] = 0;
            S_DDVCUT[i] = 0;

            C_VTAIL[i] = i_vtail*(2.*M_PI*DENS_RED[i]*NPART*DENSITY*(A[i]*LAMBDA*(RCUT_SQ + 2*RCUT*LAMBDA + 2*LAMBDA*LAMBDA)*exp(-RCUT/LAMBDA) - C[i]/3./pow(RCUT, 3) - D[i]/5./pow(RCUT, 5)));
            S_VTAIL[i] = 0;

            C_PTAIL[i] = i_vtail*(-2./3.*M_PI*DENS_RED[i]*DENSITY*DENSITY*(2.*C[i]/RCUT/RCUT/RCUT + 8.*D[i]/5./RCUT/RCUT/RCUT/RCUT/RCUT - A[i]*exp(-RCUT/LAMBDA)*(RCUT*RCUT*RCUT + 3.*RCUT*RCUT*LAMBDA + 6.*RCUT*LAMBDA*LAMBDA + 6.*LAMBDA*LAMBDA*LAMBDA)));
            S_VTAIL[i] = 0;

        } else {

            printf("\ncommon.c -> ReadInput() ERROR: Unrecognized Potential Model: POT = %c\n", POT);
            exit(EXIT_FAILURE);
        }
    }

    fscanf(fp_input, "%lf %*[^\n]\n", &LOW_TOL);

    fscanf(fp_input, "%c %*[^\n]", &dummy);
    fscanf(fp_input, "%d %*[^\n]", &IANSCREEN);
    fscanf(fp_input, "%d %*[^\n]", &IANFILE);
    fscanf(fp_input, "%d %*[^\n]", &IPS);
    fscanf(fp_input, "%d %*[^\n]", &IVMD);
    fscanf(fp_input, "%d %*[^\n]", &ISHAKE);
    fscanf(fp_input, "%d %*[^\n]", &IGOFR);
    fscanf(fp_input, "%d %*[^\n]", &ICHECK);

    if (IANSCREEN < IANFILE) {

        ANALYSIS_FLAG = 1;
        IANFILE = IANSCREEN;
        printf("\ninput.c -> ReadInput() Warning: ISCREEN < IFILE not allowed. Setting IFILE = ISCREEN\n");
    }

    fclose(fp_input);

    if (pointer_flag != 0) {
        printf("\ncommon.c -> ReadInput() ERROR: Dynamical memory allocation failure: pointer_flag = %d\n", pointer_flag);
        exit(EXIT_FAILURE);
    }

    if (access(PRECONFFN, R_OK) != -1) {

        if ((fp_input = fopen(PRECONFFN, "r+")) == NULL){
            printf("\ninput.c -> ReadInput() ERROR: Unreadable configuration file: '%s'!\n", PRECONFFN);
            exit(EXIT_FAILURE);
        }

        printf("\nReading intial configuration from file: '%s'.\n\n", PRECONFFN);

        fscanf(fp_input, "%d %lf", &tnpart, &tlbox);

        if (tnpart != NPART || fabs(tlbox-LBOX)>1e-10){

            printf("common.c -> ReadInput() ERROR: Incompatible configuration file :'%s'.\ntnpart = %d, NPART = %d\ntlbox = %.15lf, LBOX = %.15lf\nExecution aborted.\n\n", PRECONFFN, tnpart, NPART, tlbox, LBOX);
            exit(EXIT_FAILURE);
        }

        for (i=0; i<NPART; i++) {

          //fscanf(fp_input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PARTPOS_TM1[i].x, &PARTPOS_TM1[i].y, &PARTPOS_TM1[i].z, &SHELLPOS_TM1[i].x, &SHELLPOS_TM1[i].y, &SHELLPOS_TM1[i].z, &PARTPOS_T[i].x, &PARTPOS_T[i].y, &PARTPOS_T[i].z, &SHELLPOS_T[i].x, &SHELLPOS_T[i].y, &SHELLPOS_T[i].z);

          // scan checkpoint.txt for Velocity Verlet
          //fscanf(fp_input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PARTPOS_T[i].x, &PARTPOS_T[i].y, &PARTPOS_T[i].z, &PARTVEL[i].x, &PARTVEL[i].y, &PARTVEL[i].z, &SHELLPOS_T[i].x, &SHELLPOS_T[i].y, &SHELLPOS_T[i].z, &SHELLVEL[i].x, &SHELLVEL[i].y, &SHELLVEL[i].z);
          //fscanf(fp_input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PARTPOS_T[i].x, &PARTPOS_T[i].y, &PARTPOS_T[i].z, &PARTVEL[i].x, &PARTVEL[i].y, &PARTVEL[i].z, &SHELLPOS_TM1[i].x, &SHELLPOS_TM1[i].y, &SHELLPOS_TM1[i].z, &SHELLPOS_T[i].x, &SHELLPOS_T[i].y, &SHELLPOS_T[i].z);
          if (EWALD == 'F') {

              if (MODE == 'P') {

                fscanf(fp_input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PARTPOS_T[i].x, &PARTPOS_T[i].y, &PARTPOS_T[i].z, &PARTVEL[i].x, &PARTVEL[i].y, &PARTVEL[i].z, &SHELLPOS_TM1[i].x, &SHELLPOS_TM1[i].y, &SHELLPOS_TM1[i].z, &SHELLPOS_T[i].x, &SHELLPOS_T[i].y, &SHELLPOS_T[i].z);
              }

              if (MODE == 'S') {

                fscanf(fp_input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PARTPOS_T[i].x, &PARTPOS_T[i].y, &PARTPOS_T[i].z, &PARTVEL[i].x, &PARTVEL[i].y, &PARTVEL[i].z, &SHELLPOS_TM1[i].x, &SHELLPOS_TM1[i].y, &SHELLPOS_TM1[i].z, &SHELLPOS_T[i].x, &SHELLPOS_T[i].y, &SHELLPOS_T[i].z);
                SHELLVEL[i].x = PARTVEL[i].x;
                SHELLVEL[i].y = PARTVEL[i].y;
                SHELLVEL[i].z = PARTVEL[i].z;
                SHELLPOS_T[i].x = PARTPOS_T[i].x;
                SHELLPOS_T[i].y = PARTPOS_T[i].y;
                SHELLPOS_T[i].z = PARTPOS_T[i].z;
              }
          }
        }

        fclose(fp_input);

    }else{

        PRECONFIG_FLAG = 1;

        printf("\nStarting from lattice configuration.\n\n");

        InitConf();
    }

    printf("Program correctly initialized!\n\n");

    printf("#################################################################\n");

    if (MODE == 'P'){

        printf("* Polarization Mode\n");

        if (SRMODE == 'S') {

            printf("Shell Relaxation using SHAKE algorithm\n\n");

        }else if (SRMODE == 'D') {

            printf("Shell Relaxation using Steepest Descent algorithm\n\n");

        }else if (SRMODE == 'C') {

            printf("Shell Relaxation using Conjugate Gradient algorithm\n\n");

        } else {

            printf("\ncommon.c -> ReadInput() ERROR: Unrecognized Shell Relaxation Mode: SRMODE = '%c'", SRMODE);
        }

    } else if (MODE == 'S'){

        printf("* Standard Mode\n");
        printf("Automatically condensing charge on the core.\n\n");

    }else{

        printf("\ncommon.c -> ReadInput() ERROR: Unrecognized mode of simulation: MODE = '%c'\n", MODE);
        exit(EXIT_FAILURE);
    }

    if (EWALD == 'T'){

        printf("* Coulomb long-range forces treated via Ewald summation method (provisional)\n\n");

    } else if (EWALD == 'F'){

        printf("* Coulomb long-range forces cut at .5*LBOX. Second order smoothing used.\n\n");

    }else{

        printf("\ncommon.c -> ReadInput() ERROR: Unrecognized option for long-range forces: EWALD = '%c'\n", EWALD);
        exit(EXIT_FAILURE);
    }

    if (POT == 'J') {

        printf("* Using Jacucci Potential\n\n");

    } else if (POT == 'C'){

        printf("* Using Ciccotti Potential\n\n");

    }else{

        printf("\ncommon.c -> ReadInput() ERROR: Unrecognized Potential Model: POT = '%c'\n", POT);
        exit(EXIT_FAILURE);
    }

    if (THERMOSTAT == 'T') {

        printf("* Equilibration Run - NVT Run\n\n");

    } else if (THERMOSTAT == 'F'){

        printf("* Sampling Microcanonical Ensemble - NVE Run\n\n");

    }else{

        printf("\ncommon.c -> ReadInput() ERROR: Unrecognized Thermostatting Option: THERMOSTAT = '%c'\n", THERMOSTAT);
        exit(EXIT_FAILURE);
    }

    printf("* Parameters of the system: [RU] ([PhU])\n");
    printf("Number of Molecules: NMOLECULES = %d\n", NMOLECULES);
    printf("Number of Atoms per Ion: NATOMSPERMOLECULE = %d\n", NATOMSPERMOLECULE);
    printf("Number of Atomic Species: NATOMSPEC = %d\n", NATOMSPEC);
    printf("Number of Particles: NPART = %d\n", NPART);
    printf("Number of pair interactions: NINTER = %d\n", NINTER);
    printf("Length of the simulation box: LBOX = %.4e (%.4e Å)\n", LBOX, LBOX*_L_CONV);
    printf("Short range cut-off radius: RCUT = %.4e (%.4e Å)\n", RCUT, RCUT*_L_CONV);
    printf("Long range cut-off radius: LR_RCUT = %.4e (%.4e Å)\n", LR_RCUT, LR_RCUT*_L_CONV);
    printf("Long range cut-off coefficient value: LR_VCUT = %.4e (%.4e erg)\n", LR_VCUT, LR_VCUT*_E_CONV);
    printf("Density: DENSITY = %.4e (%.4e Å^(-3))\n", DENSITY, DENSITY*_DENS_CONV);
    printf("Temperature: TEMP = %.4e (%.4e K)\n", TEMP, TEMP*_TEMP_CONV);
    printf("Timestep: DT = %.4e (%.4e fs)\n", DT, DT*_TIME_CONV);
    printf("Number of timesteps: NTIMESTEPS = %d\n", NTIMESTEPS);
    printf("Total Simulation time: %.4e fs\n", NTIMESTEPS*DT*_TIME_CONV);

    printf("\n* Ions Parameters:\n");
    for (i=0; i<NATOMSPEC; i++) {

        printf(" - Atom Species: NAME[%d] = %s\n", i, NAME[i]);
        printf("Number of %s per molecule: NATOMSPERSPEC[%d] = %d\n", NAME[i], i, NATOMSPERSPEC[i]);
        printf("%s Core Mass: M[%d] = %.4e (%.4e Amu)\n", NAME[i], i, M[i], M[i]*_M_CONV);
        printf("%s Shell Mass: MU[%d] = %.4e (%.4e Amu)\n", NAME[i], i, MU[i], MU[i]*_M_CONV);
        printf("%s Core Charge: Q[%d] = %.4e (%.4e e)\n", NAME[i], i, Q[i], Q[i]*_Q_CONV);
        printf("%s Shell Charge: CHI[%d] = %.4e (%.4e e)\n", NAME[i], i, CHI[i], CHI[i]*_Q_CONV);
        printf("%s Total ionic charge: Q[%d] + CHI[%d] = %.4e (%.4e e)\n", NAME[i], i, i, Q[i] + CHI[i], (Q[i] + CHI[i])*_Q_CONV);
        printf("%s Shell-Core Coupling: K[%d] = %.4e\n", NAME[i], i, K[i]);
    }

    printf("\n* Interaction Potential Parameters and cut-off Potential Values:\n");
    for (i=0; i<NINTER; i++) {

        n1 = (int) floor(i/2.);
        n2 = (int) floor((i+1)/2.);

        printf("%s%s interaction parameters:\n", NAME[n1], NAME[n2]);
        printf("A[%d] = %.4e\n", i, A[i]);
        printf("C[%d] = %.4e\n", i, C[i]);
        printf("D[%d] = %.4e\n", i, D[i]);

        printf("\nC_VCUT[%d] = %.4e\n", i, C_VCUT[i]);
        printf("C_DVCUT[%d] = %.4e\n", i, C_DVCUT[i]);
        printf("C_DDVCUT[%d] = %.4e\n", i, C_DDVCUT[i]);

        printf("\nS_VCUT[%d] = %.4e\n", i, S_VCUT[i]);
        printf("S_DVCUT[%d] = %.4e\n", i, S_DVCUT[i]);
        printf("S_DDVCUT[%d] = %.4e\n", i, S_DDVCUT[i]);

        printf("\nC_VTAIL[%d] = %.4e\n", i, C_VTAIL[i]);

        printf("S_VTAIL[%d] = %.4e\n\n", i, S_VTAIL[i]);
    }

    printf("\nTotal Mass of the system: MTOT = %lf (%lf Amu)\n", MTOT, MTOT*_M_CONV);
    printf("SHAKE Tolerance: %.4e\n", LOW_TOL);
    printf("Output frequencies:\n");
    printf(" - Screen: ISCREEN = %d\n", IANSCREEN);
    printf(" - File: IFILE = %d\n", IANFILE);
    printf(" - PS Config: IPS = %d\n", IPS);
    printf(" - VMD: IVMD = %d\n", IVMD);
    printf(" - SHAKE: ISHAKE = %d\n", ISHAKE);
    printf(" - g(r): IGOFR = %d\n", IGOFR);
    printf(" - Checkpoint: ICHECK = %d\n", ICHECK);

    printf("#################################################################\n");
    printf("\n");
}

void FreePointers(void){

    free(NATOMSPERSPEC);
    free(PARTPOS_TM1);
    free(PARTPOS_T);
    free(PARTPOS_TP1);
    free(PARTMOM_T);
    free(PARTMOM_TP05);
    free(PARTMOM_TP1);
    free(PARTVEL);
    free(SHELLPOS_TM1);
    free(SHELLPOS_T);
    free(SHELLPOS_TP1);
    free(SHELLVEL);
    free(SHELLVEL_TP1);
    free(SHELLACC_TM1);

    if (EWALD == 'T') {

        free(CF);
        free(SF);
    }

    if (SRMODE == 'S') {

        free(DPHIDRHO_T);
        free(DPHIDVRHO_T);
        free(GAMMA);
    }

    if (SRMODE == 'C') {

        free(PHI);
        free(PHI_OLD);
        free(SEARCHDIR);
        free(RHO_OLD);
    }

    free(NAME);
    free(M);
    free(MU);
    free(Q);
    free(CHI);
    free(K);
    free(SIGMA);
    free(INDX);
    free(A);
    free(C);
    free(D);
    free(GCOUNT);
    free(GOFR);
    free(C_VCUT);
    free(S_VCUT);
    free(C_DVCUT);
    free(S_DVCUT);
    free(C_DDVCUT);
    free(S_DDVCUT);
    free(C_VTAIL);
    free(S_VTAIL);
    free(C_PTAIL);
    free(S_PTAIL);
    free(DENS_RED);
}
