//Oct. 2019 Minju Kum
//Calculates the mean Coulomb potential energy between two electrons in Helium atom
//integration uses brute force Monte-Carlo method
//infinity is approximated with lambda = 2.5
//random generator: mersenne_twister_64
//PDF: uniform distribution
//seed: number of ticks for short time duration

#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>
#include <time.h>
#include "lib.h"

using namespace std;

double integrand_bruteMC(double*);

int main(int argc, char const *argv[]) {

  long int n;
  double x[6], f;
  double sum_f = 0.;
  double sum_f_squared = 0.;
  double variance = 0.;
  double lambda = 2.5;
  double jacobi_determinant = pow(2.*lambda,6.);
  clock_t start, finish;

  //read in n
  cout << "Enter the number of Monte-Carlo samples" << endl;
  cin >> n;
  cout << n<<endl;

  //set up random number generator, with time duration seed
  //with uniform distribution in the range of (-lambda,lambda)
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  uniform_real_distribution<double> distribution(-lambda,lambda);

  start = clock();

  //measurement : n times
  for (int i = 0; i < n; i++) {
    //an array of random variables (x1,y1,z1,x2,y2,z2)
    for (int j = 0; j < 6; j++) {
      x[j] = distribution(generator);
    }
    f = integrand_bruteMC(x);
    sum_f += f;
    sum_f_squared += f*f;
  }

  //calculate mean value and variance
  sum_f = sum_f/(double(n));
  sum_f_squared = sum_f_squared/(double(n));
  variance = sum_f_squared - sum_f*sum_f;

  finish = clock();

  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << "exact result = " <<setw(10)<<setprecision(8)<< 5*atan(1)*atan(1)/16 <<endl;
  cout << "Monte-Carlo result = " <<setw(10)<<setprecision(8)<< jacobi_determinant*sum_f <<endl;
  cout << "error = " <<setw(10)<<setprecision(8)<< 5*atan(1)*atan(1)/16 - jacobi_determinant*sum_f <<endl;
  cout << "STD = " <<setw(10)<<setprecision(8)<< jacobi_determinant*sqrt(variance/double(n)) <<endl;
  cout << "CPU time: " << double(finish-start)/CLOCKS_PER_SEC <<endl;

  return 0;
} //end of main function


//define the integrand
double integrand_bruteMC(double* x){
  double alpha = 2.;
  double exponent1 = -2*alpha*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double exponent2 = -2*alpha*sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
  double denominator = sqrt(pow((x[3]-x[0]),2) + pow((x[4]-x[1]),2) + pow((x[5]-x[2]),2));
  if (denominator == 0) {
    return 0;
  }
  else return exp(exponent1+exponent2)/denominator;
}
