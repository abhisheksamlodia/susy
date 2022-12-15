#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "utils.h"

double calculate_correlator(double X[], int i){
    // returns the correlator for lattice site - i
    // X[] - 1d dimensionless field array - could be bosonic or fermionic
    // correlator = X[i]*X[0]

    return X[i]*X[0];
}

#endif