#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "utils.h"

void sparse_matrix(double phi[], double s[], double Ms[]){
    // returns nothing, updates dot product of fermion matrix M[] with vector s[]
    // phi[] - 1d dimensionless bosonic field array
    // s[] - the vector array to be used in dot product
    // Ms[] - stores the dot product value

    for (int i = 0; i < L; i++){
        Ms[i] = (((m_lat + r + (3.0 * g_lat * phi[i] * phi[i])) * s[i]) + 
                                                            (0.5*(1.0-r)*s[(i+1)%L]) - (0.5*(1.0+r)*s[(i-1+L)%L])) ;
    }

    return ;
}
#endif