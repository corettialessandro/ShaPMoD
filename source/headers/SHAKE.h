//
//  ShaPMoD.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef SHAKE_h
#define SHAKE_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "ConstTensor.h"
#include "Forces.h"

#include "output.h"

void ML_SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount);

//void SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[]);
//void VEC_SHAKE(struct point r_t[], struct point r_OLD[]);
//void ML_VEC_SHAKE(struct point r_t[], struct point r_OLD[]);

void TEST(struct point r_t[], struct point r_OLD[]);

#endif /* SHAKE_h */
