//
//  Analysis.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/21/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef Analysis_h
#define Analysis_h

#include <stdio.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"

#include "output.h"
#include "Forces.h"
#include "aux_func.h"

struct point Velocity(struct point r_tm1, struct point r_tp1);
struct point CMVelocity(struct point v[]);
struct point MultiCMVelocity(struct point v[]);

double EnerKin(struct point v[]);
double MultiEnerKin(struct point v[]);
double EnerPot_HM(struct point r[]);
double EnerPot_Jac(struct point r[], struct point rho[]);
double EnerPot_Cicc(struct point r[], struct point rho[]);
double EnerPot_WCA(struct point r[]);
double EnerPot_LJ(struct point r[]);
double EnerTot(struct point r[], struct point rho[], struct point v[]);
double Temperature(struct point v[]);
double MultiTemperature(struct point v[]);
double Pressure_Jac(struct point r[], struct point rho[]);
double Pressure_Cicc(struct point r[], struct point rho[]);
struct point Angular_Momentum(struct point r[], struct point v[]);
void GofR(struct point r[]);

void Analyse(int timestep, struct point r[], struct point rho[], struct point v[], char therm);

#endif /* Analysis_h */
