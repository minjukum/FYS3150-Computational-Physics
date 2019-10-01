//Sep 2019, Minju Kum
//Harmonic potential, 2 electrons solver
//utilize Jacobi's method

#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>
#include "jacobi.h"

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {

  cout.precision(5);
  int n; //number of steps
  int k, l;
  int count = 0;
  double h, threshold, w, interaction;
  double cpu_time, jacobicputime;
  double rho_end, rho_start; rho_start = 0;
  double Akl;
  double minimum_eigval;
  clock_t start, finish, jacobistart, jacobifinish;

  //read in endpoint rho, number of steps n, interaction
  cout << "Enter endpoint, number of steps" << endl;
  cin >> rho_end >> n;
  cout << "Set value for w" << endl;
  cin >> w;
  cout << "If there's an interaction between the electrons, enter 1. If not, enter 0" << endl;
  cin >> interaction;

  mat A = mat(n,n);
  mat V = mat(n,n); V.eye();
  vec rho = vec (n+2);
  vec pair_eigvec = vec(n);

  h = (rho_end - rho_start)/(n+1); //stepsize
  threshold = pow(10,-10); //threshold for 0

  //initialize matrix and vector
  initialize_harmonic_2el(A, rho, n, h, w, interaction, rho_end);

  start = clock();

  //if the maximum is greater than threshold
  //iterate find maximum & jacobi
  Akl = A(n-1, n-2); //set the last nondiagonal as a starting point
  while (abs(Akl)>threshold) {

    //find the maximum nondiagonal element in A
    find_maximum(A, k, l, Akl, n);

    jacobistart = clock();

    //execute jacobi method
    jacobi(A, V, n, Akl, k, l);

    jacobifinish = clock();

    count++;
  }

  finish = clock();
  cpu_time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  jacobicputime = (double(jacobifinish)-double(jacobistart))/CLOCKS_PER_SEC;

  cout << "Yes! diagonalization completed." << endl;
  cout << "# of iterations: " << count << endl;
  cout << "CPU time(s): " << cpu_time << endl;
  cout << "jacobi time(s): " <<jacobicputime*count<<endl;
  cout << "Results below: " << endl;

  //print the smallest eigenvalue and corresponding eigenvector
  vec d = A.diag();
  minimum_eigval = min(d);
  for (int i = 0; i < n; i++) {
    if (d(i) == minimum_eigval) {
      pair_eigvec = V.col(i);
    }
  }

  cout << "the smallest eigenvalue is: " <<endl<< minimum_eigval <<endl;
  cout << "corresponding eigenvector is: " <<endl<< pair_eigvec <<endl;

  cout << "Cheer up!" << endl;

  return 0;
}
