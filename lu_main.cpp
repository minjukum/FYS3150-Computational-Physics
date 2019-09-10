//10. Sep. 2019. written by Minju Kum
//this program solves differential equation d^2u(x)/dx^2 = -g(x)
//second derivative is approximated with the eqation below;
//d^2f(x)/dx^2 ~ (f(x+h)-2f(x)+f(x-h)) / (h^2).
//A(coefficient matrix) * y(vector to solve) = w(source term vector)
//in this case, the source term g is 100*exp(-10x)
//Method:
//LU decomposition
//prints vector y, minimum and maximum error,
//CPU time used only for the algorithm

#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>


using namespace std;
using namespace arma;
/*
void read_in_arguments(double &startpoint, double &endpoint, double &startval, double &endval, int &n)
void initialize(int n, double &h, double startpoint, double endpoint, vec &x, vec &f, vec &g)
void setup_tridiagonal_vectors(int n, vec &A_a, vec &A_b, vec &A_c)
void solve_tridiagonal_gaussian(int n, vec &A_b_tilde, vec A_b, vec &g_tilde, vec g, vec A_a, vec A_c, vec &f)
void get_relative_error(int n, vec f, vec r, vec &rel_error, vec &min_max_rel_error)
*/
int main(int argc, char const *argv[]) {
  int n;
  double startpoint, endpoint, startval, endval, h;
  double time_elapsed;
  clock_t start, finish;

  //read in range of x, boundary conditions, number of steps
  cout << "Enter startpoint, endpoint, start value, end value, number of steps between the two points" <<endl;
  cin >> startpoint >> endpoint >> startval >> endval >> n;

  //vectors and a matrix
  //for LU decomposition
  vec x = vec(n+2);
  mat A = mat(n,n);
  mat L = mat(n,n);
  mat U = mat(n,n);
  vec y = vec(n);
  vec z = vec(n);
  vec w = vec(n);
  double sum;

  //for error calculation
  vec r = vec(n+2);
  vec rel_error = vec(n);
  vec min_max_rel_error = vec(2);

  //define step size h
  h = (endpoint - startpoint)/(double(n+1));

  //vector x as a function of h
  for (int i = 0; i < n+2; i++) {
    x(i) = startpoint + i * h;
  }

/*
  //note:these can't be done here
  y(0) = startval;
  y(n+1) = endval;
*/

  //make coefficient matrix
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i == j)
          A(i,j) = 2.0;
        else if (abs(i-j) == 1)
          A(i,j) = -1.0;
        else
          A(i,j) = 0.0;
      }
  }

  //define w
  for (int i = 0; i < n; i++) {
    w(i) = h * h * 100 * exp(-10 * x(i+1));
  }

  start = clock();

  //LU decompose A
  lu(L,U,A);

  //Ay = LUy = Lz = w, solve Lz = w
  z(0) = w(0);
  for (int i = 1; i < n; i++) {
    sum = double(0);
    for (int j = 0; j < i; j++) {
    sum = sum + L(i,j)*z(j);
    }
    z(i) = w(i) - sum;
  }

  //solve Uy = z
  y(n-1) = z(n-1) / U(n-1,n-1);
  for (int i = n-2; i > -1; i--) {
    sum = double(0);
    for (int j = n-1; j > i ; j--) {
      sum = sum + U(i,j)*y(j);
    }
    y(i) = (z(i)-sum) / U(i,i);
  }

  finish = clock();
  time_elapsed = double(finish-start)/CLOCKS_PER_SEC;

  //to compare with the exact value
  for (int i = 0; i < n+2; i++) {
    r(i) = 1 - (1 - exp(-10)) * x(i) - exp(-10 * x(i));
  }

  //compute the relative error
  for (int i = 0; i < n; i++) {
    rel_error(i) = log10(abs((y(i)-r(i+1))/r(i+1)));
  }

  //extract the extremum value
  min_max_rel_error << min(rel_error) << max(rel_error);

  //print the result
  cout << "y is:" <<endl<< y << "min, max relative error is:" <<endl;
  cout << min_max_rel_error <<endl<< "CPU time(sec) is: " << time_elapsed <<endl;


  return 0;
}
