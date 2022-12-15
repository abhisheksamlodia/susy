#ifndef ENERGY_H
#define ENERGY_H

#include "utils.h"

double energy(double p[], double P[]){
    // returns conjugate momentum contribution for the hamiltonian
    // p[] - dimensionless conjugate momentum for 1d dimensionless bosonic field
    // P[] - dimensionless conjugate momentum for 1d dimensionless pseudo-fermionic field

    double e = 0.0;
    for (int i = 0; i < L; i++){
        e += (0.5 * ((pow(p[i],2)) + (pow(P[i],2))));
    }

    return e;
}

#endif