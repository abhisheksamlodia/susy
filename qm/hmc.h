#ifndef HMC_H
#define HMC_H

#include "utils.h"
#include "action.h"
#include "energy.h"
#include "bosonic_force.h"
#include "fermionic_force.h"
#include "cg.h"
#include "gaussian.h"
 
double update(double phi[], double Phi[]){
    // returns the absoulte value of [(Hf-Hi)/Hf]
    // phi[] - 1d dimensionless bosonic field array
    // Phi[] - 1d dimensionless pseudo-fermionic array
    // p[] - 1d conjugate momentum array for 1d bosonic field
    // P[] - 1d conjugate momentum array for 1d pseudo-fermionic field
    // Hi, Hf - initial and final hamiltonians
    // Sb, Sf - variables to store the bosonic and fermionic actions 

    double phi_old[L], Phi_old[L], p[L], P[L], s[L], f[L], F[L];
    double Hi, Hf, Sb, Sf;

    // copy current field configuration incase of rejected step
    for (int i = 0; i < L; i++){
        phi_old[i] = phi[i];
        Phi_old[i] = Phi[i];
    }

    // generate gaussian random values for conjugate momenta
    for (int j = 0; j < L; j++){
        p[j] = generate();
        P[j] = generate();
    }

    // initial hamiltonian
    Hi = action(phi, Phi, Sb, Sf) + energy(p,P);

    // leap-frog algorithm
    for (int idx = 0; idx < Ntau; idx++){
        cg(phi,Phi,s);               // findd s[] vector solution
        bosonic_force(phi,s,f);      // calculate bosonic force
        fermionic_force(F,s);        // calculate fermionic force

        // update the momenta
        for (int site = 0; site < L; site++){
            p[site] -= (f[site] * DT * 0.5);
            P[site] -= (F[site] * DT * 0.5);
        }

        // update field
        for (int i = 0; i < L; i++){
            phi[i] += (DT*p[i]);
            Phi[i] += (DT*P[i]);
        }

        cg(phi,Phi,s);               // findd s[] vector solution
        bosonic_force(phi,s,f);      // calculate bosonic force
        fermionic_force(F,s);        // calculate fermionic force

        // final update
        for (int site = 0; site < L; site++){
            p[site] -= (f[site] * DT * 0.5);
            P[site] -= (F[site] * DT * 0.5);
        }

    }
    
    // final hamiltonian
    Hf = action(phi, Phi, Sb, Sf) + energy(p,P);

    // metropolis test
    if (drand48() < exp(-(Hf-Hi))){
        // accept
        accept +=1; 
        //std::cout << "ACCEPT" << std::endl;        
    }
    else{
        //reject - copy the backup values
        //std::cout << "REJECT" << std::endl;
        for (int k = 0; k < L; k++){
            phi[k] = phi_old[k];
            Phi[k] = Phi_old[k];
        }
    }

    // return the abs((Hf-Hi)/Hi)
    return abs(((Hf-Hi)/Hi));
}

#endif