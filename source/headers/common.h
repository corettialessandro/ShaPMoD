//
//  common.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef common_h
#define common_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "parameters.h"
#include "structs.h"

#include "InitConf.h"

extern struct point *PARTPOS_TM1;
extern struct point *PARTPOS_T;
extern struct point *PARTPOS_TP1;
extern struct point *PARTVEL;
extern struct point *PARTMOM_T;
extern struct point *PARTMOM_TP05;
extern struct point *PARTMOM_TP1;
extern struct point *SHELLPOS_TM1;
extern struct point *SHELLPOS_T;
extern struct point *SHELLPOS_TP1;
extern struct point *SHELLVEL;
extern struct point *SHELLVEL_TP1;
extern struct point *SHELLACC_TM1;
extern struct point *SHELLACC_T;
extern struct point *SHELLACC_TP1;

extern struct point *CF;
extern struct point *SF;

extern struct tensor **DPHIDRHO_T;
extern struct tensor **DPHIDVRHO_T;
extern struct point *GAMMA;

extern struct point *PHI;
extern struct point *PHI_OLD;
extern struct point *SEARCHDIR;
extern struct point *RHO_OLD;

extern char MODE;
extern char SRMODE;
extern char EWALD;
extern char POT;

extern int NMOLECULES;
extern int NATOMSPEC;
extern int *NATOMSPERSPEC;
extern int NATOMSPERMOLECULE;
extern int NPART;
extern int NINTER;

extern double DT;
extern int NTIMESTEPS;

extern double DENSITY;
extern double *DENS_RED;
extern double TEMP;

extern double R_CUT;
extern double LBOX;
extern double RCUT;
extern double RCUT_SQ;
extern double LR_RCUT;
extern double ALPHA;

extern int *INDX;

extern char **NAME;
extern double *M;
extern double *MU;
extern double *Q;
extern double *CHI;
extern double *K;
extern double *SIGMA;

extern double MTOT;

extern double B;
extern double LAMBDA;

extern double *A;
extern double *C;
extern double *D;

extern double LR_VCUT;
extern double EW_VCUT;
extern double *C_VCUT;
extern double *S_VCUT;
extern double LR_DVCUT;
extern double EW_DVCUT;
extern double *C_DVCUT;
extern double *S_DVCUT;
extern double LR_DDVCUT;
extern double EW_DDVCUT;
extern double *C_DDVCUT;
extern double *S_DDVCUT;
extern double *C_VTAIL;
extern double *S_VTAIL;
extern double *C_PTAIL;
extern double *S_PTAIL;

extern int **GCOUNT;
extern int NGCOUNT;
extern double **GOFR;
extern double DG;

extern double LOW_TOL;
extern int SR_ITERS;
extern double SR_DISCR;
extern int GET_OUT;

extern char RUNNAME[];
extern char INPUTFOL[];
extern char INPUTFN[];
extern char OUTPUTFOL[];
extern char PRECONFFN[];
extern char INITCONFFN[];

extern int IANSCREEN;
extern int IANFILE;
extern int IPS;
extern int IVMD;
extern int ISHAKE;
extern int IGOFR;
extern int ICHECK;

extern int DEBUG_FLAG;
extern int VERBOSE_FLAG;
extern int OUTPUT_FLAG;
extern int ANALYSIS_FLAG;
extern int PRECONFIG_FLAG;
extern int COLLIDING_FLAG;

extern char THERMOSTAT;

extern int NANCOUNT;
extern double EPOT;
extern double ITEMP;
extern double PRESS;
extern double ETOTAVG;
extern double EKINAVG;
extern double EPOTAVG;
extern double EPOLAVG;
extern double TEMPAVG;
extern double PRESSAVG;
extern double ETOTAVGSQ;
extern double EKINAVGSQ;
extern double EPOTAVGSQ;
extern double EPOLAVGSQ;
extern double TEMPAVGSQ;
extern double PRESSAVGSQ;
extern struct point ANMOMAVG;

void ReadInput(void);
void FreePointers(void);

#endif /* common_h */
