#ifndef CG_H
#define CG_H

#include "utils.h"
#include "sparse_matrix.h"
#include "sparse_matrix_transpose.h"

# define tol 1.0e-9              // tolerance for the residue
# define MAX_ITER 100000            // maximum number of iterations

double vector_dot_product(double x[L], double y[L]){
	// returns dot product of the vector x[] and vector y[]
	double val = 0.0;
	for (int i = 0; i < L; i++){
		val += (x[i]*y[i]);
	}
	return val;
}

double vector_matrix_dot_product(double M[L][L], double x[L], int index){
	// returns dot product of a matrix with vector with a particular row of the matrix
	double val = 0.0;
	for (int i = 0; i < L; i++){
		val += (M[index][i] * x[i]);
	}
	return val;
}

double norm_vector(double x[L]){
	// returns magnitude of the vector x[]
	double val = 0.0;
	for (int i = 0; i < L; i++){
		val += (x[i]*x[i]);
	}
	return sqrt(val);
}

void cg(double phi[], double Phi[], double s[]){
	// returns nothing, finds vector s[] for the MtMs = Phi, linear algebra problem
	// phi[] - 1d dimensionless bosonic field array
    // Phi[] - 1d dimensionless pseudo-fermionic field array 

    double Ms[L], MtMs[L];                                           // check the definition in sparse_matrix.h and sparse_matrix_transpose.h

	double res[L], p_res[L];                                         // rk, pk
	double res_norm, rr_old, rr, alpha_val, palpha, beta_val;        // norm of residue vector, rkrk value, alpha_value, pk * alpha_k
	int iter = 0;                                                    // iteration counter for conjugate gradient steps
	palpha = 0.0;
	
	// initialize the initial solution array
	for (int i = 0; i < L; i++){
		s[i] = 0.0;
	}
	
    sparse_matrix(phi,s,Ms);
    sparse_matrix_transpose(phi,Ms,MtMs);

	// calculating rk and pk
	for (int i = 0; i < L; i++){
		// rk = Ax - b, pk = -rk
		res[i] = (MtMs[i] - Phi[i]);
		p_res[i] = -res[i];
	}
	
	// norm of rk
	res_norm = norm_vector(res);
	
	// conjugate algorithm loop
	while ((res_norm > tol) && (iter < MAX_ITER)){
		// alpha[k] calculation, alpha[k] = MtMs[k]
        sparse_matrix(phi,p_res,Ms);
        sparse_matrix_transpose(phi,Ms,MtMs);

		// rr = rk * rk
		rr_old = vector_dot_product(res, res);
			
		// calculate alpha_value
		palpha = vector_dot_product(p_res, MtMs);
		alpha_val = rr_old/palpha;
			
		for (int i = 0; i < L; i ++){
			// new solution
			s[i] += (alpha_val * p_res[i]);
			
			// new residue
			res[i] += (alpha_val * MtMs[i]);
		}	
		
		// calculate beta value
		rr = vector_dot_product(res, res);
		beta_val = rr/rr_old;
			
		for (int i = 0; i < L; i ++){
			// new p_residue
			p_res[i] = (-1.0 * res[i]) + (beta_val * p_res[i]);
		}
		
		// new norm
		res_norm = norm_vector(res);
			
		// update the iteration counter
		iter += 1;		
	}
    return ;
}

#endif