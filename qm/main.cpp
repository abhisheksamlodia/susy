#include "hmc.h"
#include "gaussian.h"
#include "correlator.h"
#include "utils.h"
#include "action.h"
#include "sparse_matrix.h"
#include <time.h>

using namespace std;

// Timestamp // 
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), " Date : %Y-%m-%d  Time : %X ", &tstruct);
    return buf;
}

int main(void){
    std::cout << "START" << currentDateTime() << std::endl;
    std::cout << "Lattice Mass = " << m_lat << "\t Lattice Coupling Constant = " << g_lat << std::endl ; 

    std::ofstream file1, file2, file3;
    file1.open("../data/action.txt");
    file2.open("../data/bos_corr.txt");
    file3.open("../data/ferm_corr.txt");

    double phi[L], Phi[L], s[L], Ms[L];
    double Sb,Sf,S,H;

    // initializing the field the configuration -- random gaussian numbers
    for (int i = 0; i < L; i ++){
        phi[i] = generate();
        Phi[i] = generate();
    }

    accept = 0;
    // Monte Carlo Sweeps
    srand48(12); // seed for random number generation
    
    for (int sweep = 0; sweep < TOTAL; sweep++){
        H = update(phi, Phi);
        //cout << "|dH/H| = " << H << "\t exp(-dH/H) = " << exp(-H) << endl;
        // print results after a certain sweeps
        if (sweep % GAP == 0){
            S = action(phi, Phi, Sb, Sf);
            file1 << Sb << "\t" << Sf << "\t" << S << endl; 

            cg(phi, Phi, s);
            sparse_matrix(phi, s, Ms);

            for (int i = 0; i < L; i++){
                if (i!=(L-1)){
                    file2 << phi[i]*phi[0] << "\t";
                    file3 << s[i]*Ms[0]<< "\t";
                }
                else{
                    file2 << phi[i]*phi[0] << "\n";
                    file3 << s[i]*Ms[0] << "\n";
                }
            }

            /*for (int tau = 0; tau < L; tau++){
                double sum = 0.0;
                for (int i = 0; i < L; i++){
                    for (int j = 0; j < L; j++){
                        if (abs((i-j))==tau){
                            sum += (s[j]*Ms[i])/((double)L);
                        }
                        else{
                            continue;
                        }
                    }
                }
                if (tau!=(L-1)){
                    file3 << sum << "\t";
                }
                else{
                    file3 << sum << "\n";
                }

            }*/
        }
    }

    file1.close();
    file2.close();
    file3.close();

    std::cout << "END" << currentDateTime() << std::endl;
    std::cout << "Acceptance = " << 100.0*((double)accept / (double)TOTAL) << "%" <<std::endl;

    return 0;
}