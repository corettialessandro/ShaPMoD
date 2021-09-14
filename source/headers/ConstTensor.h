//
//  ConstTensor.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef ConstTensor_h
#define ConstTensor_h

#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "structs.h"
#include "cells.h"

#include "aux_func.h"

struct tensor ConstTens_Jac(struct point rho[], struct point r[],  int i, int k);
struct tensor ConstTens_Cicc(struct point rho[], struct point r[],  int i, int k);

#endif /* ConstTensor_h */
