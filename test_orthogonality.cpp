//Sep 2019, Minju Kum
//Eigenvalue solver
//Jacobi's method - test if eigenvalues are orthogonal

#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>
#include "jacobi.h"

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {

  int n = 5; //number of steps
  int k, l;
  int count = 0;
  double h, threshold, cpu_time;
  double Akl;
  clock_t start, finish;

  mat A = mat(n,n);
  mat V = mat(n,n); V.eye();

  h = double(1)/double(n+1); //stepsize
  threshold = pow(10,-10); //threshold for 0

  //initialize matrix and vector
  initialize_beam(A, n, h);

  //if the maximum is greater than threshold
  //iterate find maximum & execute jacobi

  start = clock();

  Akl = A(n-1, n-2);
  while (abs(Akl)>threshold) {

    //find the maximum nondiagonal element in A
    find_maximum(A, k, l, Akl, n);

    //execute jacobi method
    jacobi(A, V, n, Akl, k, l);

    count++;
  }

  finish = clock();
  cpu_time = (double(finish)-double(start))/CLOCKS_PER_SEC;

  cout << "Yes! diagonalization completed." << endl;
  cout << "# of iterations: " << count << endl;
  cout << "CPU time(s): " << cpu_time << endl;
  cout << "Results below: " << endl;
  cout << "Diagonalized A is: " <<endl<< A << endl;
  cout << "Eigenvectors(col) are: " <<endl<< V << endl;

  //unit test
  cout << "let's see if the eigenvectors are orthogonal."<<endl;

  mat T = mat(n,n); T = V.t();
  mat I = mat(n,n);

  I = V * T;

  cout << "V_transposed * V is: " <<endl<< I <<endl;
  cout << "Cheer up!" << endl;

  return 0;
}
