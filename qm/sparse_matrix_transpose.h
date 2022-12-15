#ifndef SPARSE_MATRIX_TRANSPOSE_H
#define SPARSE_MATRIX_TRANSPOSE_H

#include "utils.h"

void sparse_matrix_transpose(double phi[], double Ms[], double MtMs[]){
    // returns nothing, updates dot product of transpose of fermion matrix and the vector obtained from dot
    // product of the fermion matrix and the vector
    // phi[] - 1d dimensionless bosonic field array
    // Ms[] - vector obtained from dot product of fermion matrix and the vector
    // MtMs[] - array to store the final dot product


    for (int i = 0; i < L; i ++){
        MtMs[i] = (((m_lat + r + (3.0 * g_lat * phi[i] * phi[i])) * Ms[i]) + 
                                                            (0.5*(1.0-r)*Ms[(i-1+L)%L]) - (0.5*(1.0+r)*Ms[(i+1)%L]));
    }

    return ;
}

#endif