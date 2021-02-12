//
//  output.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 1/13/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#ifndef output_h
#define output_h

#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "aux_func.h"
#include "Analysis.h"

void SetOutput(void);

void Write_MD_setup(void);

void Write_Analysis(int timestep, double ETot, double EKin, double EPot, double EPol, char therm, double Temp, double Press, struct point CMVelocity, struct point AnMom, int SHAKE_IT, double SHAKE_DIS);
void Write_GofR(int timestep, struct point r[]);

void Write_PSConfig(int timestep, struct point r[], struct point rho[], struct point v[], struct point sv[]);
void Write_InitConfig(struct point r[], struct point rho[], struct point v[]);
void Write_PartPositions(struct point r[]);
void Write_ShellPositions(struct point rho[]);
void Write_PartVelocities(struct point v[]);
void Write_ShellVelocities(struct point sv[]);
void Write_Trajectory(struct point r[], struct point rho[]);

void Write_Elapsed_timeperstep(int timestep, double elapsed_time);
void Write_Gamma(int timestep, int iteration, struct point gamma[], struct point r_tp1[], struct point rho_t[], struct point rho_tm1[]);
void Write_SHAKE_output(int timestep, int iteration, double discr, double kdiscr);
void Write_S(int timestep, int iteration, struct point rho[], struct point r[]);

void Checkpoint(int timestep, struct point r_t[], struct point rho_t[], struct point v_t[], struct point vrho_t[]);

#endif /* output_h */
