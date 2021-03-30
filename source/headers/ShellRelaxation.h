//
//  ShellRelaxation.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 4/6/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef ShellRelaxation_h
#define ShellRelaxation_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "ConstTensor.h"
#include "Forces.h"
#include "Analysis.h"
#include "aux_func.h"

#include "output.h"

void SHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], int timestep, int ccount);
void BSHAKE(struct point rho_t[], struct point rho_OLD[], struct point r_t[], struct point r_tp1[], struct point vrho_t[], struct point vrho_OLD[], int timestep, int ccount, int timeMD);
void SteepestDescent(struct point rho[], struct point r[]);
void ConjugateGradient(struct point rho[], struct point r[]);

#endif /* ShellRelaxation_h */
