#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void initialize_beam(mat&, vec&, int, double);
void initialize_harmonic_1el(mat&, vec&, int, double, double);
void initialize_harmonic_2el(mat&, vec&, int, double, double, double, double);
void find_maximum(mat, int&, int&, double&, int);
void jacobi(mat&, mat&, int, double, int, int);
void printout(int, mat, mat);

#endif
