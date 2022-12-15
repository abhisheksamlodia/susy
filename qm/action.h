#ifndef ACTION_H
#define ACTION_H

#include "utils.h"
#include "potential.h"
#include "cg.h"

double action(double phi[], double Phi[], double &Sb, double &Sf){
    // returns total action of the system
    // phi[] - 1d dimensionless bosonic field array
    // Phi[] - 1d dimensionless pseudo-fermionic field array
    // Sb, Sf - address of bosonic and fermionic part of the actions

    double pot[L], s[L];      // initialize the array for the potential and sparse matrix vector solution
    potential(phi, pot);      // find the potential array
    cg(phi,Phi,s);

    Sb =  0.0;     // initialize the action values
    Sf = 0.0;
    for (int i = 0; i < L; i++){
        // loop over all lattice sites
        Sb += ((-0.125 * phi[i] * (phi[(i+2)%L] - (2.0 * phi[i]) + phi[(i-2+L)%L])) + (0.5 * pot[i] * pot[i]));
        Sf += (0.5 * Phi[i] * s[i]);
    }

    return (Sb+Sf);
}

#endif
