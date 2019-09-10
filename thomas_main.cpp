//10. Sep. 2019. written by Minju Kum
//this program solves differential equation d^2u(x)/dx^2 = -g(x)
//second derivative is approximated with the eqation below;
//d^2f(x)/dx^2 ~ (f(x+h)-2f(x)+f(x-h)) / (h^2).
//A(coefficient matrix) * f(vector to solve) = g(source term vector)
//in this case, the source term g is 100*exp(-10x)
//Method:
//Thomas algorithm
//prints vector f, minimum and maximum error,
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
  vec x = vec(n+2);
  vec f = vec(n+2);
  vec g = vec(n+2);
  vec g_tilde = vec(n+2);
  vec A_a = vec(n-1);
  vec A_b = vec(n);
  vec A_b_tilde = vec(n);
  vec A_a_over_b_tilde = vec(n-1);
  vec A_c = vec(n-1);

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


  //apply boundary conditions
  f(0) = startval;
  f(n+1) = endval;

  //discretized function g (multiplied by h^2)
  for (int i = 0; i < n+2; i++) {
    g(i) = h * h * 100 * exp(-10 * x(i));
  }

  //diagonal and subdiagonal vectors from A
  for (int i = 0; i < n-1; i++) {
    A_a(i) = double(-1);
  }
  for (int i = 0; i < n; i++) {
    A_b(i) = double(2);
  }
  for (int i = 0; i < n-1; i++) {
    A_c(i) = double(-1);
  }


  start = clock();

  //forward substitution(special for this case)
  A_b_tilde(0) = A_b(0);
  g_tilde(1) = g(1);
  for (int i = 0; i < n-1; i++) {
    A_b_tilde(i+1) = double(i+3)/double(i+2);
    g_tilde(i+2) = g(i+2) + g_tilde(i+1)/A_b_tilde(i);
  }

  //backward substitution(special for this case)
  f(n) = g_tilde(n) / A_b_tilde(n-1);
  for (int i = n-1; i > 0; i--) {
    f(i) = (g_tilde(i)+f(i+1))/A_b_tilde(i-1);
  }

  finish = clock();
  time_elapsed = double(finish-start)/CLOCKS_PER_SEC;

  //to compare with the exact value
  for (int i = 0; i < n+2; i++) {
    r(i) = 1 - (1 - exp(-10)) * x(i) - exp(-10 * x(i));
  }

  //compute the relative error
  for (int i = 0; i < n; i++) {
    rel_error(i) = log10(abs((f(i+1)-r(i+1))/r(i+1)));
  }

  //extract the extremum value
  min_max_rel_error << min(rel_error) << max(rel_error);

  //print the result
  cout << "f is:" <<endl<< f << "min, max relative error is:" <<endl;
  cout << min_max_rel_error <<endl<< "CPU time(sec) is: " << time_elapsed <<endl;


  return 0;
}
