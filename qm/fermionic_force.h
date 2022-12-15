#ifndef FERMIONIC_FORCE_H
#define FERMIONIC_FORCE_H

#include "utils.h"

void fermionic_force(double F[], double s[]){
    // returns nothing, updates magnitude of the fermionic force that is stored in the array
    // s[] - solution of sparse marix problem using conjugate gradient method -- Look at https://arxiv.org/pdf/hep-lat/0006013

    for (int i = 0; i < L; i++){
        F[i] = s[i];
    }
    return ;
}

#endif