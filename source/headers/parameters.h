//
//  parameters.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

//Version of the program
#define _REL 0
#define _VER 5
#define _SUBVER 4

//Physical system parameters
#define _LATTICE SC     //Lattice type for initial configuration (SC)

//Temperature deviation tolerance
#define _TEMP_TOL .05

//Physical constants and conversion Factors
#define _KB 1                       //Boltzmann constant (eV/K)
#define _COULOMB 215.3209185162     //Coulomb constant (amu/e^2)

//Fixed units
#define _E_CONV 3.38e-13    //Energy units (erg)
#define _L_CONV .317        //Length units (Å)
#define _M_CONV 1.          //Mass units (Amu)
#define _Q_CONV 1.          //Charge units (e)

//Derived units
#define _TEMP_CONV 2448.1248855429      //Temperature units (K)
#define _DENS_CONV 31.3922333041        //Density units (Å^(-3))
#define _TIME_CONV 7.026280076          //Time units (fs)
#define _PRESS_CONV 1.06105748567737e12 //Pressure units (Pa)

//Shell Relaxation Parameters
#define _UP_TOL 1.e11           //Upper tolerance for shell relaxation
#define _MAX_ITER 1e4           //Maximum number of allowed iterations for reaching convergence
#define _MAX_ATT 2              //Maximum number of attempts to get out from a local minimum
#define _SR_CG_LOW_TOL 1e-11     //Lower Tolerance for shell relaxation - Conjugate Gradient
#define _SR_SHAKE_LOW_TOL 5e-11 //Lower Tolerance for shell relaxation - Shake

//Ewald Parameters
#define _ALPHA 6.

//Analysis
#define _NHIST 200

//Debug Control
#define _D_INITCONFIG 0
#define _D_FSTEP 0
#define _D_TOT_FORCES 0
#define _D_PSCONFIG 0
#define _D_ENERGY 0
#define _D_FORCES 0
#define _D_DIST 0
#define _D_CONSTR 0
#define _D_CG 0
#define _D_SHAKE 0
#define _D_TENSOR 0
#define _D_STUCK 0
#define _D_THERMOSTAT 0

//Verbose Control
#define _V_SHAKE 1

//Output Control
#define _O_SHAKE 0

#define _MAX_STR_LENGTH 100

// Magnetic field
#define B0 0
//#define INITIALIZED 1e10

// Over-relaxation parameter
#define SOR 0.33333 //0.333333

//Lennard_Jones parameters
#define LJ_SIGMA 3.405
#define LJ_EPSILON 1.

#endif /* parameters_h */
