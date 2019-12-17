/*
Nov.2019 Minju Kum
-This program simulates Sun-Mercury system in 2d with velocity Verlet algorithm
-Gm_sun is gravitational constant G multiplied by mass of the Sun,
 which is approximated as 4*pi*pi with the unit (1AU)^3/(1year)^2
-The Sun is stationary
-General relativistic correction F_g * 3l^2/(rc)^2 was added to the gravitational force F_g
 l is Mercury's angular momentum per unit mass, r is distance to the sun, and c is speed of light
-Units of position and velocity is in AU and AU/Year

   Sun
   ...    *     *                   *
.///////.  Mercury
.//^//^/.     .              *
  .///.   *          *                         *


*/

#include <iostream>
#include <cmath>
#include <iomanip>
//#include <fstream>
#include <time.h>

using namespace std;

//ofstream outfile;

void writefile(double, double, double, double, double);

int main(int argc, char const *argv[]) {

  //variables
  long int n; //number of integration points
  double t_fin,dt,t;
  int is_relativistic;
  double x,y,x1,y1,vx,vy;
  double r,l_sq,c,correction;
  double Gm_sun,Gm_mer;
  double shift;
  clock_t start, finish;

  //read in final time, number of integration points, initial conditions
  cout << "Enter the final time and number of integration points : ";
  cin >> t_fin >> n;
  cout<< "Is it relativistic? (1/0) : ";
  cin >> is_relativistic;

  //initial ctn of the Mercury
  x = 0.3075;
  y = 0.;
  vx = 0.;
  vy = 12.44;

  //set constants
  dt = t_fin/n;
  Gm_sun = 4*acos(-1.0)*acos(-1.0); // = 4*pi*pi
  c = 63241.; // speed of light in AU/Year
  l_sq = (x*vy-y*vx)*(x*vy-y*vx);

  //velocity Verlet
  //initialize
  t = 0.0;
  r = sqrt(x*x + y*y);

  start = clock();
//  outfile.open("ss-perihelion-precession.txt"); // warning: the file is very large

  //time loop
  while (t <= t_fin) {
    //save the original positions
    x1 = x;
    y1 = y;

    //calculate correction
    correction = 3.*l_sq/(r*r*c*c);

    //change the positions
    x += dt*vx + 0.5*dt*dt*(Gm_sun*(-x)/(r*r*r))*(1.+ correction*is_relativistic);
    y += dt*vy + 0.5*dt*dt*(Gm_sun*(-y)/(r*r*r))*(1.+ correction*is_relativistic);

    //we need new r and new correction
    double r_new = sqrt(x*x + y*y);
    double correction_new = 3.*l_sq/(r_new*r_new*c*c);

    //change the velocities
    vx += 0.5*dt * Gm_sun*((1.+ correction_new*is_relativistic)*(-x)/(r_new*r_new*r_new) + (1.+ correction*is_relativistic)*(-x1)/(r*r*r));
    vy += 0.5*dt * Gm_sun*((1.+ correction_new*is_relativistic)*(-y)/(r_new*r_new*r_new) + (1.+ correction*is_relativistic)*(-y1)/(r*r*r));

    r = r_new;

    t += dt;

    //write to the file
//    writefile(t, x, y, vx, vy);

    //print the shifted angle near perihelion
    if (r <= 0.3075000001) {
      shift = atan(y/x)*(180./acos(-1))*3600.; //in arcsec
      cout << setprecision(9);
      cout << t << " years passed, distance is " << r << ", " << shift << " arcsec shifted" << endl;
    }
  }
//  outfile.close();
  finish = clock();

  cout <<setprecision(4)<< "CPU time: " << double(finish-start)/CLOCKS_PER_SEC << "s" << endl;

  return 0;
}


// void writefile(double t, double x, double y, double vx, double vy)
// {
//   outfile << setiosflags(ios::showpoint);
//   outfile << setw(15) << setprecision(8) << t;
//   outfile << setw(15) << setprecision(8) << x;
//   outfile << setw(15) << setprecision(8) << y;
//   outfile << setw(15) << setprecision(8) << vx;
//   outfile << setw(15) << setprecision(8) << vy << endl;
// }
