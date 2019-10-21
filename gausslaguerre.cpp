//Oct. 2019 Minju Kum
//Calculates the mean Coulomb potential energy between two electrons in Helium atom
//integration uses Gauss-Laguerre methods

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "lib.h"
#define EPS 3.0e-14
#define MAXIT 10

using namespace std;

double integrand_lag(double,double,double,double,double,double);
void gauss_laguerre(double*,double*,int,double);
double gammln(double);

int main(int argc, char const *argv[]) {
  int n;
  double pi = acos(-1.);

  cout << "Enter number of integration points n" << endl;
  cin >> n;

  //set up mesh points and weights
  double* theta = new double [n];
  double* wle_theta = new double [n];
  double* phi = new double [n];
  double* wle_phi = new double [n];
  double* r = new double [n+1];
  double* wla_r = new double [n+1];

  gauss_legendre(0, pi, theta, wle_theta, n);
  gauss_legendre(0, 2.*pi, phi, wle_phi, n);
  gauss_laguerre(r, wla_r, n, 2.);

  //integration
  double sum_gaulag = 0.;
  for (int i = 1; i < n+1; i++) {
    for (int j = 1; j < n+1; j++) {
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          for (int m = 0; m < n; m++) {
            for (int o = 0; o < n; o++) {
              sum_gaulag += wla_r[i]*wla_r[j]*wle_theta[k]*wle_theta[l]*wle_phi[m]*wle_phi[o]*integrand_lag(r[i],theta[k],phi[m],r[j],theta[l],phi[o]);
              }
            }
          }
        }
      }
    }


  //print
  cout << "exact result: " << 5*atan(1)*atan(1)/16 << endl;
  cout << "numerical result: " << sum_gaulag << endl;
  cout << "error: " << 5*atan(1)*atan(1)/16 - sum_gaulag << endl;

  delete [] theta;
  delete [] wle_theta;
  delete [] phi;
  delete [] wle_phi;
  delete [] r;
  delete [] wla_r;

  return 0;
}



//define the integrand
double integrand_lag(double r1,double theta1,double phi1,double r2,double theta2,double phi2){

  double numerator = sin(theta1)*sin(theta2);
  double denominator_sq = r1*r1 + r2*r2 - 2.*r1*r2*(cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2));
  if (denominator_sq >= pow(10,-10)){
    return numerator/(1024.*sqrt(denominator_sq));
  }
  else {
    return 0;
  }
}


//set up mesh points and weights
void gauss_laguerre(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}


double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#undef EPS
#undef MAXIT
