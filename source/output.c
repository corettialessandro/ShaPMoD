//
//  output.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 1/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include "output.h"

void SetOutput(void){

    int i = 0;
    char answ = 'n';

    char command[_MAX_STR_LENGTH];

    sprintf(command, "mkdir -p %s", OUTPUTFOL);

    if ((i = system(command)) != 0){ //The directory cannot be created (It already exists or worse)

        sprintf(command, "touch %stemp", OUTPUTFOL);

        if ((i = system(command)) != 0) { //It is impossible to create a file into the directory (then worse!)

            printf("\noutput.c -> SetOutput() ERROR: The output directory cannot be created and it does not already exists. Unidentified problem.\nExecution aborted.\n\n");
            exit(EXIT_FAILURE);

        } else { //The directory already exists

            printf("Output folder '%s' already present. Clearing previous output files.\n", OUTPUTFOL);

            printf("Removing all the content of the folder '%s'.\nAre you sure you want to continue (y or [n])?\n", OUTPUTFOL);
            scanf("%c", &answ);

            if (answ != 'y'){

                printf("Program terminated!\n");
                exit(EXIT_SUCCESS);

            }else{

                sprintf(command, "rm -r %s*", OUTPUTFOL);

                if ((i = system(command)) != 0){ //Clearing previous putput files

                    printf("\noutput.c -> SetOutput() ERROR: Previous output files cannot be deleted.\nExecution aborted.\n\n");
                    exit(EXIT_FAILURE);
                }

                printf("Output folder '%s' correctly cleared.\n", OUTPUTFOL);

                sprintf(command, "mkdir %sSHAKE/", OUTPUTFOL);

                if ((i = system(command)) != 0){ //Creating child directory

                    printf("\noutput.c -> SetOutput() ERROR: Impossible to create a child directory.\nExecution aborted.\n\n");
                    exit(EXIT_FAILURE);
                }

                sprintf(command, "mkdir %sPSconfig/", OUTPUTFOL);

                if ((i = system(command)) != 0){ //Creating child directory

                    printf("\noutput.c -> SetOutput() ERROR: Impossible to create a child directory.\nExecution aborted.\n\n");
                    exit(EXIT_FAILURE);
                }

                OUTPUT_FLAG = 1;
                printf("Output child folder 'PSconfig/' correctly created.\n");
                printf("Output child folder 'SHAKE/' correctly created.\n");
            }
        }

    }else{

        printf("Output folder '%s' correctly created.\n", OUTPUTFOL);

        sprintf(command, "mkdir %sSHAKE/", OUTPUTFOL);

        if ((i = system(command)) != 0){ //Creating child directory

            printf("\noutput.c -> SetOutput() ERROR: Impossible to create a child directory.\nExecution aborted.\n\n");
            exit(EXIT_FAILURE);
        }

        sprintf(command, "mkdir %sPSconfig/", OUTPUTFOL);

        if ((i = system(command)) != 0){ //Creating child directory

            printf("\noutput.c -> SetOutput() ERROR: Impossible to create a child directory.\nExecution aborted.\n\n");
            exit(EXIT_FAILURE);
        }

        OUTPUT_FLAG = 2;
        printf("Output child folder 'SHAKE/' correctly created.\n");
        printf("Output child folder 'PSconfig/' correctly created.\n");
    }
}

void Write_MD_setup(void){

    int i, n1, n2;

    FILE *fp_md_setup_out;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sMD_Setup.txt", OUTPUTFOL);

    if ((fp_md_setup_out = fopen(outputpath, "w+")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_md_setup_out, "\nMD Code for Polarizable Systems (ver. %d.%d.%d).\n\n", _REL, _VER, _SUBVER);

    if (DEBUG_FLAG == 1 && VERBOSE_FLAG == 1) {

        fprintf(fp_md_setup_out, "Execution in DEBUG and VERBOSE mode: DEBUG_FLAG = %d, VERBOSE_FLAG = %d.\n", DEBUG_FLAG, VERBOSE_FLAG);

    }else if (DEBUG_FLAG == 1){

        fprintf(fp_md_setup_out, "Execution in DEBUG mode: DEBUG_FLAG = %d.\n", DEBUG_FLAG);

    }else if (VERBOSE_FLAG == 1){

        fprintf(fp_md_setup_out, "Execution in VERBOSE mode: VERBOSE_FLAG = %d.\n", VERBOSE_FLAG);

    }else{

        fprintf(fp_md_setup_out, "Standard 'silent' execution.\n");
    }

    fprintf(fp_md_setup_out, "\nReading system parameters from input file: '%s'\n", INPUTFN);

    if (OUTPUT_FLAG == 1) {

        fprintf(fp_md_setup_out, "Output folder '%s' correctly cleared.\n", OUTPUTFOL);
        fprintf(fp_md_setup_out, "Output child folder 'PSconfig/' correctly created.\n");
        fprintf(fp_md_setup_out, "Output child folder 'SHAKE/' correctly created.\n");

    }else if (OUTPUT_FLAG == 2){

        fprintf(fp_md_setup_out, "\nOutput folder '%s' correctly created.\n", OUTPUTFOL);
        fprintf(fp_md_setup_out, "Output child folder 'PSconfig/' correctly created.\n");
        fprintf(fp_md_setup_out, "Output child folder 'SHAKE/' correctly created.\n");
    }

    if (ANALYSIS_FLAG == 1) {

        fprintf(fp_md_setup_out, "\ninput.c -> ReadInput Warning: ISCREEN < IFILE not allowed. Setting IFILE = ISCREEN\n");
    }

    if (PRECONFIG_FLAG == 0) {

        fprintf(fp_md_setup_out, "\nReading intial configuration from file: '%s'.\n\n", PRECONFFN);

    }else if (PRECONFIG_FLAG == 1){

        fprintf(fp_md_setup_out, "\nStarting from lattice configuration.\n\n");
    }

    fprintf(fp_md_setup_out, "Program correctly initialized!\n\n");

    fprintf(fp_md_setup_out, "#################################################################\n");

    if (MODE == 'P'){

        fprintf(fp_md_setup_out, "* Polarization Mode\n");

        if (SRMODE == 'S') {

            fprintf(fp_md_setup_out, "Shell Relaxation using SHAKE algorithm\n\n");

        }else if (SRMODE == 'D') {

            fprintf(fp_md_setup_out, "Shell Relaxation using Steepest Descent algorithm\n\n");

        }else if (SRMODE == 'C') {

            fprintf(fp_md_setup_out, "Shell Relaxation using Conjugate Gradient algorithm\n\n");

        }

    } else if (MODE == 'S'){

        fprintf(fp_md_setup_out, "* Standard Mode\n");
        fprintf(fp_md_setup_out, "Automatically condensing charge on the core.\n\n");
    }

    if (EWALD == 'T'){

        fprintf(fp_md_setup_out, "* Coulomb long-range forces treated via Ewald summation method (provisional)\n\n");

    } else if (EWALD == 'F'){

        fprintf(fp_md_setup_out, "* Coulomb long-range forces cut at .5*LBOX. Second order smoothing used.\n\n");
    }

    if (POT == 'J') {

        fprintf(fp_md_setup_out, "* Using Jacucci Potential\n\n");

    } else if (POT == 'C'){

        fprintf(fp_md_setup_out, "* Using Ciccotti Potential\n\n");
    }

    if (THERMOSTAT == 'T') {

        fprintf(fp_md_setup_out, "* Equilibration Run - NVT Run\n\n");

    } else if (THERMOSTAT == 'F'){

        fprintf(fp_md_setup_out, "* Sampling Microcanonical Ensemble - NVE Run\n\n");
    }

    fprintf(fp_md_setup_out, "* Parameters of the system: [RU] ([PhU])\n");
    fprintf(fp_md_setup_out, "Number of Molecules: NMOLECULES = %d\n", NMOLECULES);
    fprintf(fp_md_setup_out, "Number of Atoms per Ion: NATOMSPERMOLECULE = %d\n", NATOMSPERMOLECULE);
    fprintf(fp_md_setup_out, "Number of Atomic Species: NATOMSPEC = %d\n", NATOMSPEC);
    fprintf(fp_md_setup_out, "Number of Particles: NPART = %d\n", NPART);
    fprintf(fp_md_setup_out, "Number of pair interactions: NINTER = %d\n", NINTER);
    fprintf(fp_md_setup_out, "Length of the simulation box: LBOX = %.4e (%.4e Å)\n", LBOX, LBOX*_L_CONV);
    fprintf(fp_md_setup_out, "Short range cut-off radius: RCUT = %.4e (%.4e Å)\n", RCUT, RCUT*_L_CONV);
    fprintf(fp_md_setup_out, "Long range cut-off radius: LR_RCUT = %.4e (%.4e Å)\n", LR_RCUT, LR_RCUT*_L_CONV);
    fprintf(fp_md_setup_out, "Long range cut-off coefficient value: LR_VCUT = %.4e (%.4e erg)\n", LR_VCUT, LR_VCUT*_E_CONV);
    fprintf(fp_md_setup_out, "Density: DENSITY = %.4e (%.4e Å^(-3))\n", DENSITY, DENSITY*_DENS_CONV);
    fprintf(fp_md_setup_out, "Temperature: TEMP = %.4e (%.4e K)\n", TEMP, TEMP*_TEMP_CONV);
    fprintf(fp_md_setup_out, "Timestep: DT = %.4e (%.4e fs)\n", DT, DT*_TIME_CONV);
    fprintf(fp_md_setup_out, "Number of timesteps: NTIMESTEPS = %d\n", NTIMESTEPS);
    fprintf(fp_md_setup_out, "Total Simulation time: %.4e fs\n", NTIMESTEPS*DT*_TIME_CONV);

    fprintf(fp_md_setup_out, "\n* Ions Parameters:\n");
    for (i=0; i<NATOMSPEC; i++) {

        fprintf(fp_md_setup_out, " - Atom Species: NAME[%d] = %s\n", i, NAME[i]);
        fprintf(fp_md_setup_out, "Number of %s per molecule: NATOMSPERSPEC[%d] = %d\n", NAME[i], i, NATOMSPERSPEC[i]);
        fprintf(fp_md_setup_out, "%s Core Mass: M[%d] = %.4e (%.4e Amu)\n", NAME[i], i, M[i], M[i]*_M_CONV);
        fprintf(fp_md_setup_out, "%s Shell Mass: MU[%d] = %.4e (%.4e Amu)\n", NAME[i], i, MU[i], MU[i]*_M_CONV);
        fprintf(fp_md_setup_out, "%s Core Charge: Q[%d] = %.4e (%.4e e)\n", NAME[i], i, Q[i], Q[i]*_Q_CONV);
        fprintf(fp_md_setup_out, "%s Shell Charge: CHI[%d] = %.4e (%.4e e)\n", NAME[i], i, CHI[i], CHI[i]*_Q_CONV);
        fprintf(fp_md_setup_out, "%s Total ionic charge: Q[%d] + CHI[%d] = %.4e (%.4e e)\n", NAME[i], i, i, Q[i] + CHI[i], (Q[i] + CHI[i])*_Q_CONV);
        fprintf(fp_md_setup_out, "%s Shell-Core Coupling: K[%d] = %.4e\n", NAME[i], i, K[i]);
    }

    fprintf(fp_md_setup_out, "\n* Interaction Potential Parameters:\n");
    for (i=0; i<NINTER; i++) {

        n1 = (int) floor(i/2.);
        n2 = (int) floor((i+1)/2.);

        fprintf(fp_md_setup_out, "%s%s interaction parameters:\n", NAME[n1], NAME[n2]);
        fprintf(fp_md_setup_out, "A[%d] = %.4e\n", i, A[i]);
        fprintf(fp_md_setup_out, "C[%d] = %.4e\n", i, C[i]);
        fprintf(fp_md_setup_out, "D[%d] = %.4e\n", i, D[i]);

        fprintf(fp_md_setup_out, "\nC_VCUT[%d] = %.4e\n", i, C_VCUT[i]);
        fprintf(fp_md_setup_out, "C_DVCUT[%d] = %.4e\n", i, C_DVCUT[i]);
        fprintf(fp_md_setup_out, "C_DDVCUT[%d] = %.4e\n", i, C_DDVCUT[i]);

        fprintf(fp_md_setup_out, "\nS_VCUT[%d] = %.4e\n", i, S_VCUT[i]);
        fprintf(fp_md_setup_out, "S_DVCUT[%d] = %.4e\n", i, S_DVCUT[i]);
        fprintf(fp_md_setup_out, "S_DDVCUT[%d] = %.4e\n", i, S_DDVCUT[i]);

        fprintf(fp_md_setup_out, "\nC_VTAIL[%d] = %.4e\n", i, C_VTAIL[i]);

        fprintf(fp_md_setup_out, "S_VTAIL[%d] = %.4e\n\n", i, S_VTAIL[i]);
    }

    fprintf(fp_md_setup_out, "\nTotal mass of the system: MTOT = %lf (%lf Amu)\n", MTOT, MTOT*_M_CONV);
    fprintf(fp_md_setup_out, "SHAKE Tolerance: %.4e\n", LOW_TOL);
    fprintf(fp_md_setup_out, "Output frequencies:\n");
    fprintf(fp_md_setup_out, " - Screen: ISCREEN = %d\n", IANSCREEN);
    fprintf(fp_md_setup_out, " - File: IFILE = %d\n", IANFILE);
    fprintf(fp_md_setup_out, " - PS Config: IPS = %d\n", IPS);
    fprintf(fp_md_setup_out, " - VMD: IVMD = %d\n", IVMD);
    fprintf(fp_md_setup_out, " - SHAKE: ISHAKE = %d\n", ISHAKE);
    fprintf(fp_md_setup_out, " - g(r): IGOFR = %d\n", IGOFR);
    fprintf(fp_md_setup_out, " - Checkpoint: ICHECK = %d\n", ICHECK);

    fclose(fp_md_setup_out);
}

void Write_Analysis(int timestep, double ETot, double EKin, double EPot, double EPol, char therm, double Temp, double Press, struct point CMVelocity, struct point AnMom, int SHAKE_IT, double SHAKE_DIS){

    FILE *fp_analysis_out;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sAnalysis.txt", OUTPUTFOL);

    if ((fp_analysis_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Analysis_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_analysis_out, "%d\t%.15e\t%.10e\t%.10e\t%.10e\t%c%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%d\t%.10e\n", timestep, ETot, EKin, EPot, 100*EPol/EPot, therm, Temp, Press, CMVelocity.x, CMVelocity.y, CMVelocity.z, AnMom.x, AnMom.y, AnMom.z, SHAKE_IT, SHAKE_DIS);
    fflush(fp_analysis_out);
    fclose(fp_analysis_out);
}

void Write_GofR(int timestep, struct point r[]){

    FILE *fp_gofr_out;

    int i, h, n1, n2;
    double R, Vbin, d3h, Nig;

    char outputpath[_MAX_STR_LENGTH];

    GofR(r);

    if (timestep % (5*IGOFR) == 0 && NGCOUNT != 0) {

        for (i=0; i<NINTER; i++) {

            n1 = (int) floor(i/2.);
            n2 = (int) floor((i+1.)/2.);

            sprintf(outputpath, "%sGofR_%s%s.txt", OUTPUTFOL, NAME[n1], NAME[n2]);

            if ((fp_gofr_out = fopen(outputpath, "w+")) == NULL){

                printf("\noutput.c -> Write_GofR() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
                exit(EXIT_FAILURE);
            }

            for (h=0; h<_NHIST; h++) {

                R = DG*(h+.5);
                d3h = 3*h*h + 3*h + 1;
                Vbin = d3h*DG*DG*DG;
                Nig = 4./3.*M_PI*Vbin*DENS_RED[i]*DENSITY;
                GOFR[i][h] = (double)GCOUNT[i][h]/(NGCOUNT*(NPART-1.)*Nig);

                fprintf(fp_gofr_out, "%lf\t%lf\n", R*_L_CONV, GOFR[i][h]);
            }

            fclose(fp_gofr_out);
        }
    }
}

void Write_PSConfig(int timestep, struct point r[], struct point rho[], struct point v[], struct point sv[]){

    int y, z, tlog;
    if (timestep == 0) {

        Write_InitConfig(r, rho, v);

    } else {

        if (MODE == 'S') {

//            for (z=0;z<=7;z++) {
//                for (y=1;y<=9;y++) {
//                    tlog = y*(int)(pow(10.,z));
//                    if (timestep==0 || timestep==tlog) {
                        if ((timestep%20) == 0){
                          Write_PartPositions(r,timestep);
                          Write_PartVelocities(v,timestep);
                        }
//                    }
//                }
//            }

        } else if (MODE == 'P'){

            Write_PartPositions(r,timestep);
            Write_ShellPositions(rho);
            Write_PartVelocities(v,timestep);
            Write_ShellVelocities(sv);
        }
    }
}

void Write_InitConfig(struct point r[], struct point rho[], struct point v[]){

    FILE *fp_psconfig_out;

    int i;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sPSconfig/Initial_Configuration.txt", OUTPUTFOL);

    if ((fp_psconfig_out = fopen(outputpath, "w+")) == NULL){

        printf("\noutput.c -> Write_InitConfig() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_psconfig_out, "#Initial Configuration\n");
    fprintf(fp_psconfig_out, "#\n#part\tpart_indx\tcoremass\tshellmass\tcorecharge\tshellcharge\ttotalcharge\tcoreposxyz\tcorevelxyz\tshellposxyz\n");

    for (i=0; i<NPART; i++) {

        fprintf(fp_psconfig_out, "%d\t%d\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", i, INDX[i], M[INDX[i]], MU[INDX[i]], Q[INDX[i]], CHI[INDX[i]], Q[INDX[i]]+CHI[INDX[i]], r[i].x, r[i].y, r[i].z, v[i].x, v[i].y, v[i].z, rho[i].x, rho[i].y, rho[i].z);
    }

    fclose(fp_psconfig_out);
}

void Write_PartPositions(struct point r[], int time){

    FILE *fp_partpos_out0;
    FILE *fp_partpos_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/r_neg_t%d.txt", OUTPUTFOL,time);
    sprintf(outputpath1, "%sPSconfig/r_pos_t%d.txt", OUTPUTFOL,time);

    if ((fp_partpos_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_PartPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_partpos_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_PartPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_partpos_out0, "%.10e\t%.10e\t%.10e\n", r[i].x, r[i].y, r[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_partpos_out1, "%.10e\t%.10e\t%.10e\n", r[i].x, r[i].y, r[i].z);
        }
    }

    fclose(fp_partpos_out0);
    fclose(fp_partpos_out1);
}

void Write_LogPartPositions(struct point r[], int timestep){

    FILE *fp_partpos_out0;
    FILE *fp_partpos_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/r_neg_logt%d.txt", OUTPUTFOL,timestep);
    sprintf(outputpath1, "%sPSconfig/r_pos_logt%d.txt", OUTPUTFOL,timestep);

    if ((fp_partpos_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_PartPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_partpos_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_PartPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_partpos_out0, "%.10e\t%.10e\t%.10e\n", r[i].x, r[i].y, r[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_partpos_out1, "%.10e\t%.10e\t%.10e\n", r[i].x, r[i].y, r[i].z);
        }
    }

    fclose(fp_partpos_out0);
    fclose(fp_partpos_out1);
}

void Write_ShellPositions(struct point rho[]){

    FILE *fp_shellpos_out0;
    FILE *fp_shellpos_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/rho_neg.txt", OUTPUTFOL);
    sprintf(outputpath1, "%sPSconfig/rho_pos.txt", OUTPUTFOL);

    if ((fp_shellpos_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_ShellPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_shellpos_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_ShellPositions() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_shellpos_out0, "%.10e\t%.10e\t%.10e\n", rho[i].x, rho[i].y, rho[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_shellpos_out1, "%.10e\t%.10e\t%.10e\n", rho[i].x, rho[i].y, rho[i].z);
        }
    }

    fclose(fp_shellpos_out0);
    fclose(fp_shellpos_out1);
}

void Write_PartVelocities(struct point v[], int time){

    FILE *fp_partvel_out0;
    FILE *fp_partvel_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/v_neg_t%d.txt", OUTPUTFOL,time);
    sprintf(outputpath1, "%sPSconfig/v_pos_t%d.txt", OUTPUTFOL,time);

    if ((fp_partvel_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_PartVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_partvel_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_PartVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_partvel_out0, "%.10e\t%.10e\t%.10e\n", v[i].x, v[i].y, v[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_partvel_out1, "%.10e\t%.10e\t%.10e\n", v[i].x, v[i].y, v[i].z);
        }
    }

    fclose(fp_partvel_out0);
    fclose(fp_partvel_out1);
}

void Write_LogPartVelocities(struct point v[], int timestep){

    FILE *fp_partvel_out0;
    FILE *fp_partvel_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/v_neg_logt%d.txt", OUTPUTFOL,timestep);
    sprintf(outputpath1, "%sPSconfig/v_pos_logt%d.txt", OUTPUTFOL,timestep);

    if ((fp_partvel_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_PartVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_partvel_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_PartVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_partvel_out0, "%.10e\t%.10e\t%.10e\n", v[i].x, v[i].y, v[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_partvel_out1, "%.10e\t%.10e\t%.10e\n", v[i].x, v[i].y, v[i].z);
        }
    }

    fclose(fp_partvel_out0);
    fclose(fp_partvel_out1);
}

void Write_ShellVelocities(struct point sv[]){

    FILE *fp_shellvel_out0;
    FILE *fp_shellvel_out1;

    int i;

    char outputpath0[_MAX_STR_LENGTH];
    char outputpath1[_MAX_STR_LENGTH];

    sprintf(outputpath0, "%sPSconfig/sv_neg.txt", OUTPUTFOL);
    sprintf(outputpath1, "%sPSconfig/sv_pos.txt", OUTPUTFOL);

    if ((fp_shellvel_out0 = fopen(outputpath0, "a")) == NULL){

        printf("\noutput.c -> Write_ShellVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath0);
        exit(EXIT_FAILURE);
    }

    if ((fp_shellvel_out1 = fopen(outputpath1, "a")) == NULL){

        printf("\noutput.c -> Write_ShellVelocities() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath1);
        exit(EXIT_FAILURE);
    }

    for (i=0; i<NPART; i++) {

        if (INDX[i] == 0) {

            fprintf(fp_shellvel_out0, "%.10e\t%.10e\t%.10e\n", sv[i].x, sv[i].y, sv[i].z);

        } else if (INDX[i] == 1) {

            fprintf(fp_shellvel_out1, "%.10e\t%.10e\t%.10e\n", sv[i].x, sv[i].y, sv[i].z);
        }
    }

    fclose(fp_shellvel_out0);
    fclose(fp_shellvel_out1);
}

void Write_Trajectory(struct point r[], struct point rho[]){

    FILE * fp_vmd_out;

    int i, indx_i;
    char core_name[10], shell_name[10];
    struct point zero = {0};

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sTrajectory.xyz", OUTPUTFOL);

    if ((fp_vmd_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Write_VMD() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_vmd_out, "%d\n", 2*NPART);
    fprintf(fp_vmd_out, "\n");

    for (i=0; i<NPART; i++) {

        indx_i = INDX[i];

        sprintf(core_name, "%s", NAME[indx_i]);
        sprintf(shell_name, "Sh_%s", NAME[indx_i]);

        fprintf(fp_vmd_out, "%s\t%lf\t%lf\t%lf\n%s\t%lf\t%lf\t%lf\n", core_name, d_rirj(r[i], zero).x*_L_CONV, d_rirj(r[i], zero).y*_L_CONV, d_rirj(r[i], zero).z*_L_CONV, shell_name, d_rhoirj(rho[i], r[i], zero).x*_L_CONV, d_rhoirj(rho[i], r[i], zero).y*_L_CONV, d_rhoirj(rho[i], r[i], zero).z*_L_CONV);
//        fprintf(fp_vmd_out, "%s\t%lf\t%lf\t%lf\n%s\t%lf\t%lf\t%lf\n", core_name, Distance(r[i], zero).x*_L_CONV, Distance(r[i], zero).y*_L_CONV, Distance(r[i], zero).z*_L_CONV, shell_name, Distance(rho[i], zero).x*_L_CONV, Distance(rho[i], zero).y*_L_CONV, Distance(rho[i], zero).z*_L_CONV);
    }

    fclose(fp_vmd_out);
}

void Write_Elapsed_timeperstep(int timestep, double elapsed_time){

    FILE *fp_time;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sTime.txt", OUTPUTFOL);

    if ((fp_time = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Write_Elapsed_timeperstep() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_time, "%d\t%lf\n", timestep, elapsed_time);

    fclose(fp_time);
}

void Write_Gamma(int timestep, int iteration, struct point gamma[], struct point r_tp1[], struct point rho_t[], struct point rho_tm1[]){

    FILE *fp_gamma_out;

    int k;
    double kk;
    struct point gamma_theor, delta_gamma;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sSHAKE/gamma_out.txt", OUTPUTFOL);

    if ((fp_gamma_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Write_Gamma() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    if (timestep == 1 && iteration == 0) fprintf(fp_gamma_out, "#timestep\titeration\tconstr_index\tg_num\tg_theor\t∆gamma\n");

    for (k=0; k<NPART; k++) {

        kk = K[INDX[k]];

        gamma_theor.x=(r_tp1[k].x + rho_tm1[k].x - 2*rho_t[k].x)/kk;
        gamma_theor.y=(r_tp1[k].y + rho_tm1[k].y - 2*rho_t[k].y)/kk;
        gamma_theor.z=(r_tp1[k].z + rho_tm1[k].z - 2*rho_t[k].z)/kk;

        delta_gamma.x = gamma[k].x-gamma_theor.x;
        delta_gamma.y = gamma[k].y-gamma_theor.y;
        delta_gamma.z = gamma[k].z-gamma_theor.z;

        fprintf(fp_gamma_out, "%d\t%d\t%d\t%.7e\t%.7e\t%.7e\n", timestep, iteration, 3*k+0, gamma[k].x, gamma_theor.x, delta_gamma.x);
        fprintf(fp_gamma_out, "%d\t%d\t%d\t%.7e\t%.7e\t%.7e\n", timestep, iteration, 3*k+1, gamma[k].y, gamma_theor.y, delta_gamma.y);
        fprintf(fp_gamma_out, "%d\t%d\t%d\t%.7e\t%.7e\t%.7e\n", timestep, iteration, 3*k+2, gamma[k].z, gamma_theor.z, delta_gamma.z);
    }

    fclose(fp_gamma_out);
}

void Write_SHAKE_output(int timestep, int iteration, double discr, double kdiscr){

    FILE *fp_shake_out;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sSHAKE/SHAKE_out_%.04d.txt", OUTPUTFOL, timestep);

    if ((fp_shake_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Write_SHAKE_output() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    fprintf(fp_shake_out, "%d\t%.25e\t%.1lf\n", iteration, discr, kdiscr);

    fclose(fp_shake_out);
}

void Write_S(int timestep, int iteration, struct point rho[], struct point r[]) {

    FILE *fp_s_out;

    int i;

    char outputpath[_MAX_STR_LENGTH];

    sprintf(outputpath, "%sSHAKE/S_out.txt", OUTPUTFOL);

    if ((fp_s_out = fopen(outputpath, "a")) == NULL){

        printf("\noutput.c -> Write_S() ERROR: File %s not found.\nExecution aborted.\n\n", outputpath);
        exit(EXIT_FAILURE);
    }

    if (timestep == 1 && iteration == 0) fprintf(fp_s_out, "#timestep\titeration\tparticle\ts_x\ts_y\ts_z\n");

    for (i=0; i<NPART; i++) {

        fprintf(fp_s_out, "%d\t%d\t%d\t%.10e\t%.10e\t%.10e\n", timestep, iteration, i, rho[i].x - r[i].x, rho[i].y - r[i].y, rho[i].z - r[i].z);
    }

    fclose(fp_s_out);
}

void Checkpoint(int timestep, struct point r_t[], struct point rho_t[], struct point rho_tm1[], struct point v_t[], struct point vrho_t[]){

    FILE *fp_checkpoint;

    int i;

//    struct point zero = {0,0,0};
    struct point tr_t, trho_t, trho_tm1, tv_t, tvrho_t;

    if (timestep == 0) {

        if ((fp_checkpoint = fopen(INITCONFFN, "w+")) == NULL){

            printf("\noutput.c -> Checkpoint() ERROR: Impossible to create file %s.\nExecution aborted.\n\n", INITCONFFN);
            exit(EXIT_FAILURE);
        }

    } else {

        if ((fp_checkpoint = fopen(PRECONFFN, "w+")) == NULL){

            printf("\noutput.c -> Checkpoint() ERROR: Impossible to create file %s.\nExecution aborted.\n\n", PRECONFFN);
            exit(EXIT_FAILURE);
        }
    }

    fprintf(fp_checkpoint, "%d\t%.15lf\n", NPART, LBOX);

    for (i=0; i<NPART; i++) {

// //        tr_tm1 = d_rirj(r_tm1[i], zero);
// //        tr_tm1 = Distance(r_tm1[i], zero);
//         tr_tm1 = r_tm1[i];
// //        trho_tm1 = d_rhoirj(rho_tm1[i], r_tm1[i], zero);
// //        trho_tm1 = Distance(rho_tm1[i], zero);
//         trho_tm1 = rho_tm1[i];
// //        tr_t = d_rirj(r_t[i], zero);
// //        tr_t = Distance(r_t[i], zero);
//         tr_t = r_t[i];
// //        trho_t = d_rhoirj(rho_t[i], r_t[i], zero);
// //        tr_t = Distance(rho_t[i], zero);
//         trho_t = rho_t[i];
//
//         fprintf(fp_checkpoint, "%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n", tr_tm1.x, tr_tm1.y, tr_tm1.z, trho_tm1.x, trho_tm1.y, trho_tm1.z, tr_t.x, tr_t.y, tr_t.z, trho_t.x, trho_t.y, trho_t.z);

      // create temporary variables for particle&shell position&velocity
      tr_t = r_t[i];
      //tr_tm1 = r_t[i];
      tv_t = v_t[i];
      trho_t = rho_t[i];
      trho_tm1 = rho_tm1[i];
      tvrho_t = vrho_t[i];
      // create save.txt
      fprintf(fp_checkpoint, "%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n", tr_t.x, tr_t.y, tr_t.z, trho_t.x, trho_t.y, trho_t.z, trho_tm1.x, trho_tm1.y, trho_tm1.z, tv_t.x, tv_t.y, tv_t.z, tvrho_t.x, tvrho_t.y, tvrho_t.z);
    }

    fprintf(fp_checkpoint, "\nAverage Thermodynamics Variables of the simulation:\n");
    fprintf(fp_checkpoint, "<E> = %.4e\t∆E/<E> = %.4e\t<K> = %.4e\t<U> = %.4e\t100<Up>/<U> = %.4e\t<T> = %.4e\t<P> = %.4e\n", ETOTAVG, sqrt(fabs(ETOTAVGSQ - ETOTAVG*ETOTAVG))/fabs(ETOTAVG), EKINAVG, EPOTAVG, 100*EPOLAVG/EPOTAVG, TEMPAVG, PRESSAVG);

    fclose(fp_checkpoint);

}
