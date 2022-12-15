#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "utils.h"

void potential(double phi[], double pot[]){
    // returns nothing, updates the potential array
    // phi[] - 1d dimensionless bosonic field array
    // pot[] - dimensionless potential array
    // p_pow - power of 1d bosonic field for the potential
    
    for (int i = 0; i < L; i ++){
        pot[i] = (((m_lat + r) * phi[i]) - (0.5*r*(phi[(i+1)%L] + phi[(i-1+L)%L])) + (g_lat * phi[i] * phi[i] * phi[i])); 
    }

    return ;
}

#endif