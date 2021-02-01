//
//  FirstStep.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef FirstStep_h
#define FirstStep_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "Analysis.h"
#include "Forces.h"
#include "SHAKE.h"
#include "ShellRelaxation.h"
#include "output.h"

void FirstStep_St(void);
void FirsrStep_St_Ew(void);

void FirstStep_Pol(void);
void FirstStep_Pol_Ew(void);

#endif /* FirstStep_h */
