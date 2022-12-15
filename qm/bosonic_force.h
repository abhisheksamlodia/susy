#ifndef BOSONIC_FORCE_H
#define BOSONIC_FORCE_H

# include "utils.h"
# include "potential.h"
# include "sparse_matrix.h"

void bosonic_force(double phi[], double s[], double f[]){
    // returns nothing, updates bosonic force to the array f[]
    // phi[] - 1d dimensionless bosonic field array
    // s[] - solution vector for sparse matrix problem to calculate forces

    double pot[L], Ms[L];                  // initializing potential and M*s vector, M is fermion matrix
    potential(phi, pot);                   // find the potential
    sparse_matrix(phi, s, Ms);             // find the M*s vector for given s vector
    
    for (int i = 0; i < L; i++){
        f[i] = (((-0.25 * (phi[(i+2)%L] - (2.0 * phi[i]) + phi[(i-2+L)%L]))) + 
                    (((m_lat + r + ((double)p_pow * g_lat * pow(phi[i], (p_pow-1))))*pot[i]) - 
                        (0.5*r*pot[(i+1)%L]) - (0.5*r*pot[(i-1+L)%L])) +
                        ((double)p_pow * (double) (p_pow-1) * g_lat * pow(phi[i], (p_pow-2)) * s[i] * Ms[i]));
    }

    return ;
}

#endif