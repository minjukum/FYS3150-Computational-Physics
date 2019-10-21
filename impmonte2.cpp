//Oct. 2019 Minju Kum
//Calculates the mean Coulomb potential energy between two electrons in Helium atom
//integration uses Monte-Carlo method with importance sampling
//random generator: mersenne_twister_64
//PDF: exponential & uniform distribution
//seed: number of ticks of short time duration

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <time.h>
#include "lib.h"

using namespace std;

double integrand_impMC(double*);

int main(int argc, char const *argv[]) {

  int n;
  double pi = acos(-1.);
  double x[6], f;
  double sum_f = 0.;
  double sum_f_squared = 0.;
  double variance = 0.;
  double jacobi_determinant = 4*pow(acos(-1.),4.);
  clock_t start, finish;

  cout << "Enter the number of Monte-Carlo samples" << endl;
  cin >> n;

  //set up random number generator with time duration seed
  //set up distribution for r with exponential PDF of mean = 1/4
  //and uniform PDF for theta, phi
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  exponential_distribution<double> r_distribution (4.0);
  uniform_real_distribution<double> theta_distribution(0.,pi);
  uniform_real_distribution<double> phi_distribution(0.,2.*pi);

  start = clock();

  //sample n times
  for (int i = 0; i < n; i++) {
    //arrange r1, theta1, phi1 and r2, theta2, phi2
    for (int j = 0; j < 2; j++) {
      x[0+3*j] = r_distribution(generator);     // 0 < r1, r2 < inf
      x[1+3*j] = theta_distribution(generator); // 0 < theta1, theta2 < pi
      x[2+3*j] = phi_distribution(generator);   // 0 < phi1, phi2 < 2pi
    }
    f = integrand_impMC(x);
    sum_f += f;
    sum_f_squared += f*f;
  }

  sum_f = sum_f/(double(n));
  sum_f_squared = sum_f_squared/(double(n));
  variance = sum_f_squared - sum_f*sum_f;

  finish = clock();

  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << "exact result = " <<setw(10)<<setprecision(8)<< 5*atan(1)*atan(1)/16 <<endl;
  cout << "Monte-Carlo result = " <<setw(10)<<setprecision(8)<< jacobi_determinant*sum_f <<endl;
  cout << "error = " <<setw(10)<<setprecision(8)<< 5*atan(1)*atan(1)/16 - jacobi_determinant*sum_f <<endl;
  cout << "STD = " <<setw(10)<<setprecision(8)<< jacobi_determinant*sqrt(variance/double(n)) <<endl;
  cout << "CPU time: " << double(finish-start)/CLOCKS_PER_SEC << endl;

  return 0;
} //end of main function


//define the integrand
double integrand_impMC(double* x){

  double numerator_MC = x[0]*x[0]*x[3]*x[3]*sin(x[1])*sin(x[4]);
  double denominator_sq_MC = x[0]*x[0] + x[3]*x[3] -2.*x[0]*x[3]*(cos(x[1])*cos(x[4]) + sin(x[1])*sin(x[4])*cos(x[2]-x[5]));

  if (denominator_sq_MC >= pow(10,-20)) {
    return numerator_MC/(16.*sqrt(denominator_sq_MC));
  }
  else return 0;
}
