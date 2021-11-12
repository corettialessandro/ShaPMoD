//
//  aux_func.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef aux_func_h
#define aux_func_h

#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

struct point Distance(struct point A, struct point B);
struct point d_rirj(struct point ri, struct point rj);
struct point d_rirhoj(struct point ri, struct point rhoj, struct point rj);
struct point d_rhoirj(struct point rhoi, struct point ri, struct point rj);
struct point d_rhoirhoj(struct point rhoi, struct point ri, struct point rhoj, struct point rj);

double Inner(struct point A, struct point B);
struct point Vector(struct point A, struct point B);

double mod(struct point A);
double modsq(struct point A);

struct point RNDM_Rotate(struct point r);

double Variance(struct point v[], int n);
double OTRAvg(int N, double Estimator, double Avg);

double Thermostat(struct point r_tm1[], struct point r_t[], double temp);

void LinearConjugateGradient(double **matrix, struct point *vector_b, struct point *vector_x, int ndimension, char component);

#endif /* aux_func_h */
