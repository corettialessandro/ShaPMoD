//
//  lattice.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef lattice_h
#define lattice_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

struct point SC(int part_indx, int npart);
struct point BCC(int part_indx, int npart);

#endif /* lattice_h */
