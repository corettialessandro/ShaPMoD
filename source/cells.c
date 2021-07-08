#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include "parameters.h"
#include "structs.h"
#include "common.h"
#include "cells.h"

#define ID3D(i,j,k) (k + (j + i*N_CELLS_Y)*N_CELLS_Z);

/*** TODO ***
 * 
 * we need to define the following variables (global?)
 * N_CELLS_X : (int)(floor(LX/RCUT))
 * N_CELLS_Y : (int)(floor(LY/RCUT))
 * N_CELLS_Z : (int)(floor(LZ/RCUT))
 * CELL_SIZE_X : LX/N_CELLS_X
 * CELL_SIZE_Y : LY/N_CELLS_Y
 * CELL_SIZE_Z : LZ/N_CELLS_Z
 * 
 * we need to define the following arrays (global?)
 * int HEADS[N_CELLS_X * N_CELLS_Y * N_CELLS_Z]
 * int HOST[N_POINTS]
 * int LINKS[N_POINTS]
 * !!! These arrays need to be initialised to -1 !!!
 *
************/

int Find_Cell ( const struct point p )
{
	int i = (int)((p.x+0.5*LBOX)/CELL_SIZE_X);
	int j = (int)((p.y+0.5*LBOX)/CELL_SIZE_Y);
	int k = (int)((p.z+0.5*LBOX)/CELL_SIZE_Z);
	return ID3D(i,j,k);
}

void Find_Ind ( const int cell_label, int *i, int *j, int *k ) //CHECK THIS
{
	*i = (int)(cell_label/(N_CELLS_Y*N_CELLS_Z));
	*j = (int)(cell_label/N_CELLS_Z) - (int)(cell_label/(N_CELLS_Y*N_CELLS_Z))*N_CELLS_Y;
	*k = cell_label % N_CELLS_Z;
//	printf("%d -> %d\n",l,ID3D(*i,*j,*k));
}

void Add_Point_To_Cell ( const struct point p, const int label)
{
	int n = Find_Cell(p);
	//printf("Adding particle %d [%lf,%lf,%lf] -> [%lf,%lf,%lf] into cell %d [%d,%d,%d]\n",label,p.x,p.y,p.z,p.x+0.5*LBOX,p.y+0.5*LBOX,p.z+0.5*LBOX,n,(int)((p.x+0.5*LBOX)/CELL_SIZE_X),(int)((p.y+0.5*LBOX)/CELL_SIZE_Y),(int)((p.z+0.5*LBOX)/CELL_SIZE_Z));
	//printf("LINKS[%d] = HEADS[%d]\n",label,n);
	LINKS[label] = HEADS[n];
	//printf("HEADS[%d] = %d\n",n,label);
	HEADS[n] = label;
	//printf("HOST[%d] = %d\n",label,n);
	HOST[label] = n;
}

void Rem_Point_From_Cell (const int label)
{
	int n = HEADS[HOST[label]];
	if (n==label) HEADS[HOST[label]] = LINKS[n];
	else { while(LINKS[n]!=label) n = LINKS[n]; }
	LINKS[n] = LINKS[label];
	HOST[label] = -1;
}

void List_Of_Neighs ( const int label, int *list, int nshells )
{
// Writes neighbors' labels on list
// The first element is the number of neighbors
// and the remaining elements are the labels of the neighboring particles.
// The list also contains the particle itself!!!
	int i,j,k,ii,jj,kk,pi,pj,pk;
	int c,l;
	int n=0;
	Find_Ind(HOST[label],&i,&j,&k);
	for (ii=-nshells;ii<=nshells;ii++) {
		for (jj=-nshells;jj<=nshells;jj++) {
			for (kk=-nshells;kk<=nshells;kk++) {
				pi = (i + ii + N_CELLS_X) % N_CELLS_X;
				pj = (j + jj + N_CELLS_Y) % N_CELLS_Y;
				pk = (k + kk + N_CELLS_Z) % N_CELLS_Z;
				c = ID3D(pi,pj,pk);
				l = HEADS[c];
				while (l>=0) {
					list[n+1] = l;
					l = LINKS[l];
					n++;
				}
			}
		}
	}
	*list = n;
}
