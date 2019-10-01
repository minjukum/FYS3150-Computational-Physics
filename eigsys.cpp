//Sep 2019, Minju Kum
//buckling beam solver
//utilize armadillo function eig_sym()

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
  vec eigval;
  mat eigvec;

  h = double(1)/double(n+1); //stepsize

  //initialize matrix and vector
  double h_squared = h*h;
  for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        if (i == j)
          A(i,j) = 2/h_squared;
        else if (abs(i-j) == 1)
          A(i,j) = -1/h_squared;
        else
          A(i,j) = 0.0;
      }
  }

  start = clock();

  eig_sym(eigval, eigvec, A);

  finish = clock();
  cpu_time = (double(finish)-double(start))/CLOCKS_PER_SEC;

  cout << "Yes! diagonalization completed." << endl;
  cout << "CPU time(s): " << cpu_time << endl;
  cout << "Cheer up!" << endl;
  printout(n, eigval, eigvec);

  return 0;
}
