//
//  Block_MD.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef Block_MD_h
#define Block_MD_h

#include <stdio.h>
#include <time.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "FirstStep.h"
#include "Analysis.h"
#include "Forces.h"
#include "SHAKE.h"
#include "ShellRelaxation.h"
#include "aux_func.h"
#include "output.h"

void Block_MD_St(void);
void Block_MD_St_Ew(void);

void Block_MD_Pol(void);
void Block_MD_Pol_Ew(void);

void Block_MD_MultiMaze(void);

#endif /* Block_MD_h */
