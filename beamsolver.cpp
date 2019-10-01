//Sep 2019, Minju Kum
//buckling beam solver
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
  double h, threshold, cpu_time;
  double Akl;
  clock_t start, finish;

  //read in number of steps n
  cout << "Enter number of steps" <<endl;
  cin >> n;

  mat A = mat(n,n);
  mat V = mat(n,n); V.eye();
  vec rho = vec(n+2);

  h = double(1)/double(n+1); //stepsize
  threshold = pow(10,-10); //threshold for 0

  //initialize matrix and vector
  initialize_beam(A, rho, n, h);

  start = clock();

  //if the maximum is greater than threshold
  //iterate find maximum & jacobi
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

  //if (n>=10) print the 3 smallest eigenvalues and corresponding eigenvectors
  //else print all
  printout(n, A, V);

  cout << "Cheer up!" << endl;

  return 0;
}
