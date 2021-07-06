#ifndef cells_h
#define cells_h

#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include "parameters.h"
#include "structs.h"

void Init_Cell ( int cut );
int Find_Cell ( const struct point p );
void Find_Ind ( const int cell_label, int *i, int *j, int *k ); //CHECK THIS
void Add_Point_To_Cell ( const struct point p, const int label );
void Rem_Point_From_Cell (const int label );
void List_Of_Neighs ( const int label, int *list, int nshells );

#endif /* cells_h */