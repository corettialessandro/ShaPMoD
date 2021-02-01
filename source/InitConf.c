//
//  InitConf.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#include "InitConf.h"

void InitConf(void){
    
    int i, indx_i;
    struct point vcm = {0}, rcm = {0};
    struct point CF = {0}, SF = {0};
    double mi, alpha, v2 = 0, g = 3.*(NPART - 1.);
    
    srand48(time(0));
    
    for (i=0; i<NPART; i++) {
        
        mi = M[INDX[i]];
        
        SHELLPOS_T[i] = PARTPOS_T[i] = _LATTICE(i, NPART);

        if (COLLIDING_FLAG == 1) {
            
            PARTVEL[i].x = pow(-1, i);
            PARTVEL[i].y = pow(-1, i);
            PARTVEL[i].z = pow(-1, i);
        
        }else{
        
            PARTVEL[i].x = (double)lrand48()/(RAND_MAX+1.) - .5;
            PARTVEL[i].y = (double)lrand48()/(RAND_MAX+1.) - .5;
            PARTVEL[i].z = (double)lrand48()/(RAND_MAX+1.) - .5;
        }

        rcm.x += mi*PARTPOS_T[i].x;
        rcm.y += mi*PARTPOS_T[i].y;
        rcm.z += mi*PARTPOS_T[i].z;
        
        vcm.x += mi*PARTVEL[i].x;
        vcm.y += mi*PARTVEL[i].y;
        vcm.z += mi*PARTVEL[i].z;
        
        v2 += mi*modsq(PARTVEL[i]);
    }
    
    rcm.x /= MTOT;
    rcm.y /= MTOT;
    rcm.z /= MTOT;

    vcm.x /= MTOT;
    vcm.y /= MTOT;
    vcm.z /= MTOT;
    
    v2 /= g;
    
    alpha = sqrt(TEMP/v2);
    
    for (i=0; i<NPART; i++) {
        
        SHELLPOS_T[i].x = PARTPOS_T[i].x -= rcm.x;
        SHELLPOS_T[i].y = PARTPOS_T[i].y -= rcm.y;
        SHELLPOS_T[i].z = PARTPOS_T[i].z -= rcm.z;

        PARTVEL[i].x = (PARTVEL[i].x - vcm.x)*alpha;
        PARTVEL[i].y = (PARTVEL[i].y - vcm.y)*alpha;
        PARTVEL[i].z = (PARTVEL[i].z - vcm.z)*alpha;
    }
    
    if (DEBUG_FLAG && _D_INITCONFIG) {
        
        if (POT == 'J') {
            
            for (i=0; i<NPART; i++) {
                
                CF = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
                SF = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);
                
                printf("F[%d] = (%.4e, %.4e, %.4e)\tPhi[%d] = (%.4e, %.4e, %.4e)\n", i, CF.x, CF.y, CF.z, i, SF.x, SF.y, SF.z);
            }
            
        }else if (POT == 'C'){
            
            for (i=0; i<NPART; i++) {
                
                CF = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
                SF = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
                
                printf("F[%d] = (%.4e, %.4e, %.4e)\tPhi[%d] = (%.4e, %.4e, %.4e)\n", i, CF.x, CF.y, CF.z, i, SF.x, SF.y, SF.z);
            }
        }
    }
    
    if (DEBUG_FLAG && _D_INITCONFIG) {
        
        if (POT == 'J') {
            
            for (i=0; i<NPART; i++) {
                
                CF = CoreForce_Jac(PARTPOS_T, SHELLPOS_T, i);
                SF = ShellForce_Jac(SHELLPOS_T, PARTPOS_T, i);
                
                printf("F[%d] = (%.4e, %.4e, %.4e)\tPhi[%d] = (%.4e, %.4e, %.4e)\n", i, CF.x, CF.y, CF.z, i, SF.x, SF.y, SF.z);
            }
            
        }else if (POT == 'C'){
            
            for (i=0; i<NPART; i++) {
                
                CF = CoreForce_Cicc(PARTPOS_T, SHELLPOS_T, i);
                SF = ShellForce_Cicc(SHELLPOS_T, PARTPOS_T, i);
                
                printf("F[%d] = (%.4e, %.4e, %.4e)\tPhi[%d] = (%.4e, %.4e, %.4e)\n", i, CF.x, CF.y, CF.z, i, SF.x, SF.y, SF.z);
            }
        }
    }

    vcm.x = vcm.y = vcm.z = 0;
    v2 = 0;
    
    for (i=0; i<NPART; i++) {
        
        mi = M[INDX[i]];
        
        vcm.x += mi*PARTVEL[i].x;
        vcm.y += mi*PARTVEL[i].y;
        vcm.z += mi*PARTVEL[i].z;

        v2 += mi*modsq(PARTVEL[i]);
    }
    
    if (DEBUG_FLAG && _D_INITCONFIG) {
        
        printf("\nINITIAL CONFIGURATION:\n");
        printf("i\tindx\tM\tMu\tQ\tChi\tx\ty\tz\tvx\tvy\tvz\n");
        
        for (i=0; i<NPART; i++) {
            
            indx_i = INDX[i];
            
            printf("%d\t%d\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", i, indx_i, M[indx_i], MU[indx_i], Q[indx_i], CHI[indx_i], PARTPOS_T[i].x, PARTPOS_T[i].y, PARTPOS_T[i].z, PARTVEL[i].x, PARTVEL[i].y, PARTVEL[i].z);
        }
    }
    
    Write_PSConfig(0, PARTPOS_T, SHELLPOS_T, PARTVEL, PARTVEL);
        
    printf("Initial configuration correctly generated.\n");
    printf("Initial total momentum: vcm = (%.4e, %.4e, %.4e)\n", vcm.x/MTOT, vcm.y/MTOT, vcm.z/MTOT);
    printf("System initial temperature: T = %.4e K\n\n", v2/g*_TEMP_CONV);
}
