#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

const int L = 16;                        // number of lattice sites
const double DT = 0.1;                  // leapfrog time step
const int TOTAL = 5000000;                // number of hmc sweeps
const int GAP = 1000;                    // gap between each measurement
const double r = 1.0;                    // parameter of wilson mass operator 
const int p_pow = 3;                     // power of the bosonic field in super-potential 
const double m_phy = 10.0;               // physical mass in wilson mass operator
const double g_phy = 10.0;              // coupling constant in super-potential

int Ntau = (int)(1.0/DT);                // number of leapfrog updates
double a = (1.0/(double)L);              // lattice spacing 
double m_lat = m_phy * a;                // dimensionless mass
double g_lat = g_phy * a * a;            // dimensionless coupling constant
int accept;                              // acceptance counter	

#endif