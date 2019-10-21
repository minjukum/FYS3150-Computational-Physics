//Oct. 2019 Minju Kum
//Calculates the mean Coulomb potential energy between two electrons in Helium atom
//integration uses Gauss-Legendre methods
//gauss_legendre() is gauleg() from lib.h

#include <iostream>
#include <cmath>
#include "lib.h"

using namespace std;

double integrand(double,double,double,double,double,double);

int main(int argc, char const *argv[]) {
  int N;
  double a, b;

  cout << "Enter integration limits a b and number of integration points N" <<endl;
  cin >> a >> b >> N;

  //mesh points and weights
  double* xle = new double [N];
  double* wle = new double [N];

  //return the abcissas in x[0:N-1] and w[0:N-1]
  gauss_legendre(a, b, xle, wle, N);

  //brute force legendre
  double sum_gauleg = 0.;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        for (int l = 0; l < N; l++) {
          for (int m = 0; m < N; m++) {
            for (int n = 0; n < N; n++) {
              sum_gauleg+=wle[i]*wle[j]*wle[k]*wle[l]*wle[m]*wle[n]*integrand(xle[i],xle[j],xle[k],xle[l],xle[m],xle[n]);
            }
          }
        }
      }
    }
  }
/*
  //for plot
  double* r = new double [N];
  for (int i = 0; i < N; i++) {
    r[i] = i*b/(N-1);
  }
  double* integ = new double [N];
  for (int i = 0; i < N; i++) {
    double alpha = 2.;
    integ[i] = exp(-alpha*r[i]);
  }

  //prints
  for (int i = 0; i < N; i++) {
    cout << r[i] << endl;
  }
  cout << endl;
  for (int i = 0; i < N; i++) {
    cout << integ[i] << endl;
  }
*/
  cout << "exact result: " << 5*atan(1)*atan(1)/16 <<endl;
  cout << "numerical result: " << sum_gauleg << endl;
  cout << "error: " << 5*atan(1)*atan(1)/16 - sum_gauleg << endl;

  delete [] xle;
  delete [] wle;
//  delete [] r;
//  delete [] integ;

  return 0;
}

//define the function being integrated
double integrand(double x1,double y1,double z1,double x2,double y2,double z2){
  double alpha = 2.;
  double exponent1 = -2*alpha*sqrt(x1*x1 + y1*y1 + z1*z1);
  double exponent2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
  double denominator = sqrt(pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2));
  if (denominator == 0) {
    return 0;
  }
  else return exp(exponent1+exponent2)/denominator;
}
