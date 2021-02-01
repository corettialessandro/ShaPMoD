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
