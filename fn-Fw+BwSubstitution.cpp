//forward and backward substitution
#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {
  int n;
  double s, e, h;

  //read in number of steps
  cout << "Enter startpoint, endpoint, number of steps" <<endl;
  cin >> s >> e >> n;

  //vectors and a matrix
  vec x = vec(n-1);
  vec f = vec(n-1);
  vec g = vec(n-1);
  vec g_tilda = vec(n-1);
  vec r = vec(n-1);
  mat A = mat(n-1,n-1);
  vec A_a = vec(n-2);
  vec A_b = vec(n-1);
  vec A_b_tilda = vec(n-1);
  vec A_c = vec(n-2);

  //define step size h
  h = (e - s)/(double)n;

  //vector x as a function of h
  for (int i = 0; i < n-1; i++) {
    x(i) = s + (i+1) * h;
  }

  //vectorized function g (multiplied by h^2)
  for (int i = 0; i < n-1; i++) {
    g(i) = h * h * 100 * exp(-10 * x(i));
  }

  //make coefficient matrix
  for (int i = 0; i < n-1; i++) {
      for (int j = 0; j < n-1; j++) {
        if (i==j)
          A(i,j) = 2.0;
        else if (abs(i-j)==1)
          A(i,j) = -1.0;
        else
          A(i,j) = 0.0;
      }
  }

  //diagonal and subdiagonal vectors from A
  for (int i = 0; i < n-2; i++) {
    A_a(i) = A(i+1,i);
  }
  for (int i = 0; i < n-1; i++) {
    A_b(i) = A(i,i);
  }
  for (int i = 0; i < n-2; i++) {
    A_c(i) = A(i,i+1);
  }

  //forward substitution
  for (int i = 0; i < n-2; i++) {
    A_b_tilda(0) = A_b(0);
    g_tilda(0) = g(0);
    A_b_tilda(i+1) = A_b(i+1) - (A_a(i)*A_c(i))/A_b_tilda(i);
    g_tilda(i+1) = g(i+1) - (A_a(i)*g_tilda(i))/A_b_tilda(i);
  }
  cout << A <<endl<< A_a <<endl<< A_b <<endl<< A_b_tilda <<endl
  << A_c <<endl<< g <<endl<< g_tilda <<endl;

  //backward substitution
  for (int i = n-2; i > 0; i--) {
    f(n-2) = g_tilda(n-2)/A_b_tilda(n-2);
    f(i-1) = (g_tilda(i-1) - A_c(i-1)*f(i))/A_b_tilda(i-1);
  }

  //to compare with the exact value
  for (int i = 0; i < n-1; i++) {
    r(i) = 1 - (1 - exp(-10)) * x(i) - exp(-10 * x(i));
  }

  cout << f <<endl<< r <<endl;



  return 0;
}
