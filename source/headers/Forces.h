//
//  Forces.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef Forces_h
#define Forces_h

#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "aux_func.h"
#include "ConstTensor.h"

struct point CoreForce_Jac(struct point r[], struct point rho[], int i);
struct point ShellForce_Jac(struct point rho[], struct point r[], int i);

struct point CoreForce_Cicc(struct point r[], struct point rho[], int i);
struct point ShellForce_Cicc(struct point rho[], struct point r[], int i);

double CoreVirial_Jac(struct point r[], struct point rho[], int i);
double ShellVirial_Jac(struct point rho[], struct point r[], int i);

double CoreVirial_Cicc(struct point r[], struct point rho[], int i);
double ShellVirial_Cicc(struct point rho[], struct point r[], int i);

void Forces(struct point r[], struct point rho[]);

#endif /* Forces_h */
