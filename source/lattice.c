//
//  lattice.c
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

//For binary mixtures only

#include "lattice.h"

struct point SC(int part_indx, int npart){
    
    struct point PartPos;
    
    int i,j,k;
    double jj,l;
    int PartPerSide = (int)round(pow(npart, 1./3.));
    double LattStep = LBOX/(double)PartPerSide;
    
    if ((l=pow(PartPerSide, 3)) != npart) {
        
        printf("\nlattice.c -> SC() ERROR: Bad number of particles\n");
        printf("PartPerSide = %d\tPartPerSide**3 = %.0lf\tNpart = %d\n\n", PartPerSide, l, npart);
        exit(EXIT_FAILURE);
    }
    
    i = part_indx%PartPerSide;
    jj = round((double)part_indx/(double)PartPerSide);
    j = (int)(jj+1)%PartPerSide;
    k = (int)(jj/(double)PartPerSide)%PartPerSide;
    
    //Switching first atom
    i = (i+j)%PartPerSide;
    j = (j+k)%PartPerSide;
    
    PartPos.x = i*LattStep;
    PartPos.y = j*LattStep;
    PartPos.z = k*LattStep;
    
    return PartPos;
}

struct point BCC(int part_indx, int npart){
    
    struct point PartPos;
    
    double l;
    int PartPerSide = (int)round(pow(.5*npart, 1./3.));
    double LattStep = LBOX/(double)PartPerSide;
    
    if ((l=2*pow(PartPerSide, 3)) != npart) {
        
        printf("\nlattice.c -> BCC() ERROR: Bad number of particles\n");
        printf("PartPerSide = %d\t2*(PartPerSide**3) = %.0lf\tNpart = %d\n\n", PartPerSide, l, npart);
        exit(EXIT_FAILURE);
    }
    
    if (PartPerSide%2 == 0){
        
        printf("\nlattice.c -> BCC() ERROR: Bad number of particles\n");
        printf("Insert an odd number of molecules!\n\n");
        exit(EXIT_FAILURE);
    }

    if (INDX[part_indx] == 0){

        PartPos = SC(part_indx, (int)round(.5*npart));

    }else{

        PartPos = SC(part_indx, (int)round(.5*npart));

        PartPos.x += LattStep*.5;
        PartPos.y += LattStep*.5;
        PartPos.z += LattStep*.5;
    }

    return PartPos;
}
