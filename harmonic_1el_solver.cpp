//Sep 2019, Minju Kum
//Harmonic potential, 1 electron solver
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
  double h, threshold, cpu_time, jacobicputime;
  double rho_end, rho_start; rho_start = 0;
  double Akl;
  clock_t start, finish, jacobistart, jacobifinish;

  //read in endpoint rho, number of steps n
  cout << "Enter endpoint, number of steps" << endl;
  cin >> rho_end >> n;

  mat A = mat(n,n);
  mat V = mat(n,n); V.eye();
  vec rho = vec (n+2);

  h = (rho_end - rho_start)/(n+1); //stepsize
  threshold = pow(10,-10); //threshold for 0

  //initialize matrix and vector
  initialize_harmonic_1el(A, rho, n, h, rho_end);

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

  //if (n>10) print the 3 smallest eigenvalues and corresponding eigenvectors
  //else print all
  printout(n, A, V);

  cout << "Cheer up!" << endl;

  return 0;
}
