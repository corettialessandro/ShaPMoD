//
//  main.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

//ON GOING:
// - Steepest descent method for energy minimization
// - Conjugate gradient method for energy minimization
// - Efficiency (Force, Pressure, g(r))
// - Ewald summation for electromagnetic interactions
//
//TO DO:
// - Lattice option from input file
// - Generalize interaction indexing
// - Implementation with Velocity Verlet Algorithm
// - Implementation of modified VVA for Magnetic Field
// - Nosé-Hoover Thermostat for Magnetic Field

#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"
#include "output.h"

#include "Block_MD.h"

int main(int argc, char *argv[]) {
    
    int c;
    
    printf("\nSHAKE POLARIZATION MOLECULAR DYNAMICS (ShaPMoD)\n");
    printf("MD Code for Polarizable Systems (ver. %d.%d.%d).\n\n", _REL, _VER, _SUBVER);

    while ((c = getopt(argc, argv, "cdgi:o:r:v")) != -1) {
        switch (c) {
            case 'c':
                COLLIDING_FLAG = 1;
                break;
            case 'd':
                DEBUG_FLAG = 1;
                break;
            case 'g':
                GET_OUT = 1;
                break;
            case 'i':
                sprintf(INPUTFN, "%s", optarg);
                break;
            case 'o':
                sprintf(OUTPUTFOL, "output/%s/", optarg);
                break;
            case 'r':
                sprintf(RUNNAME, "%s", optarg);
                sprintf(INPUTFOL, "input/%s/", RUNNAME);
                sprintf(OUTPUTFOL, "output/%s/", RUNNAME);
                sprintf(RUNNAME, "%s_", optarg);
                sprintf(INPUTFN, "%s%sinput.txt", INPUTFOL, RUNNAME);
                sprintf(PRECONFFN, "%s%scheckpoint.txt", INPUTFOL, RUNNAME);
                sprintf(INITCONFFN, "%s%sinit_checkpoint.txt", INPUTFOL, RUNNAME);
                break;
            case 'v':
                VERBOSE_FLAG = 1;
                break;
            case '?':
                if (optopt == 'i') {
                    printf("Option -%c requires the name of the input file.\nExecution aborted.\n\n", optopt);
                }else if (optopt == 'o') {
                    printf("Option -%c requires the name of the output file.\nExecution aborted.\n\n", optopt);
                }else if (optopt == 'r') {
                    printf("Option -%c requires the name of the run.\nExecution aborted.\n\n", optopt);
                }else if (isprint(optopt)) {
                    printf("Unknown option '-%c'.\nExecution aborted.\n\n", optopt);
                } else {
                    printf("Unknown option character '\\x%x'.Execution aborted.\n\n", optopt);
                }
                exit(EXIT_FAILURE);
            default:
                abort();
        }
    }
    
    if (access(INPUTFN, F_OK) == -1) {
        
        printf("\nmain.c -> main() ERROR: Input file '%s' not found.\nExecution aborted.\n\n", INPUTFN);
        exit(EXIT_FAILURE);
    }
    
    if (DEBUG_FLAG == 1 && VERBOSE_FLAG == 1) {
        
        printf("Execution in DEBUG and VERBOSE mode: DEBUG_FLAG = %d, VERBOSE_FLAG = %d.\n", DEBUG_FLAG, VERBOSE_FLAG);
    
    }else if (DEBUG_FLAG == 1){
        
        printf("Execution in DEBUG mode: DEBUG_FLAG = %d.\n", DEBUG_FLAG);
    
    }else if (VERBOSE_FLAG == 1){
        
        printf("Execution in VERBOSE mode: VERBOSE_FLAG = %d.\n", VERBOSE_FLAG);

    }else{
        
        printf("Standard 'silent' execution.\n\n");
    }
    
    clock_t t_start, t_end;
    t_start = clock();
    
    SetOutput();

    printf("\nReading system parameters from input file: '%s'\n", INPUTFN);

    ReadInput();

    Write_MD_setup();
    
    printf("\n#################################################################\n");
    printf("Begin of dynamics:\n");
    printf("t\tETot\tEKin\tEPot\t100*EPol/EPot\tTemperature\tPressure\tVCMx\tVCMy\tVCMz\tSR_ITs\tSR_DISCR\n");
    
    if (EWALD == 'T') {
        
        if (MODE == 'P') {
            
            Block_MD_Pol_Ew();
            
        }else if (MODE == 'S'){
            
            Block_MD_St_Ew();
        }
        
    } else if (EWALD == 'F') {
        
        if (MODE == 'P') {
            
            Block_MD_Pol();
            
        }else if (MODE == 'S'){
            
            Block_MD_St();
        }
    }
    
    printf("#################################################################\n");
    
    printf("\nAverage Thermodynamics Variables of the simulation:\n");
    printf("<E> = %.4e\t∆E/<E> = %.4e\t<K> = %.4e\t<U> = %.4e\t100*<Up>/<U> = %.4e\t<T> = %.4e\t<P> = %.4e\n", ETOTAVG, sqrt(fabs(ETOTAVGSQ - ETOTAVG*ETOTAVG))/fabs(ETOTAVG), EKINAVG, EPOTAVG, 100*EPOLAVG/EPOTAVG, TEMPAVG, PRESSAVG);

    FreePointers();
    
    t_end = clock();
    
    printf("\nTotal execution time: %.2f (s)\n\n", (double)(t_end - t_start)/CLOCKS_PER_SEC);
    return 0;
}
