/*
Nov.2019 Minju Kum
This program simulates Sun-Earth system in 3d
with forward Euler and velocity Verlet.
Gm_sun is gravitational constant G multiplied by mass of the Sun,
which is approximated as 4*pi*pi with the unit (1AU)^3/(1year)^2.
Energies and angular momentum are divided by mass of the Earth.

   Sun
   ...    *     *         *         *
.///////.     *     Earth
.//^//^/.             .     *
  .///.   *       *                         *
   *

*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;

ofstream outfile;

void writefile(double, double, double, double, double, double, double, double, double, double, double, double, double);

int main(int argc, char const *argv[]) {

  //variables
  long int n;
  double t_fin,dt,t;
  double r,x,y,z,vx,vy,vz,ax,ay,az;
  double x1,y1,z1,r_new;
  double Gm_sun;
  double E_kinetic, E_potential, E_total;
  double angular_momentum_x, angular_momentum_y,angular_momentum_z;
  clock_t eulstart, eulfinish, verstart, verfinish;

  //read in final time, number of integration points, initial conditions
  cout << "Enter the final time and number of integration points : ";
  cin >> t_fin >> n;
//  cout << "Enter the initial position x, y, z : ";
//  cin >> x >> y >> z;
//  cout << "Enter the initial velocity v_x, v_y, v_z : ";
//  cin >> vx >> vy >> vz;

  x = 1.;
  y = 0.;
  z = 0.;
  vx = 0.;
  vy = 2.*acos(-1.);
  vz = 0.;

  //set constants
  dt = t_fin/n;
  Gm_sun = 4*acos(-1.0)*acos(-1.0); // = 4*pi*pi
  //save for later
  double x0 = x;
  double y0 = y;
  double z0 = z;
  double vx0 = vx;
  double vy0 = vy;
  double vz0 = vz;

  //forward Euler no VERLET
  //initialize
  t = 0.0;
  r = sqrt(x*x + y*y + z*z);
  E_kinetic = 0.;
  E_potential = 0.;
  E_total = 0.;
  angular_momentum_x = 0.;
  angular_momentum_y = 0.;
  angular_momentum_z = 0.;

  verstart = clock();
  outfile.open("sun-earth-verlet.txt");

  //time loop
  while (t <= t_fin) {
    //save the original positions to calculate velocities
    x1 = x;
    y1 = y;
    z1 = z;
    //change the positions
    x += dt*vx - 0.5*dt*dt*Gm_sun*x/(r*r*r);
    y += dt*vy - 0.5*dt*dt*Gm_sun*y/(r*r*r);
    z += dt*vz - 0.5*dt*dt*Gm_sun*z/(r*r*r);
    //and we also need new r
    r_new = sqrt(x*x + y*y + z*z);
    //keep changing
    vx -= 0.5*dt*Gm_sun*(x/(r_new*r_new*r_new) + x1/(r*r*r));
    vy -= 0.5*dt*Gm_sun*(y/(r_new*r_new*r_new) + y1/(r*r*r));
    vz -= 0.5*dt*Gm_sun*(z/(r_new*r_new*r_new) + z1/(r*r*r));
    r = r_new;
    t += dt;
    E_kinetic = 0.5*(vx*vx + vy*vy + vz*vz);
    E_potential = -Gm_sun/r;
    E_total = E_kinetic + E_potential;
    angular_momentum_x = y*vz-z*vy;
    angular_momentum_y = z*vx-x*vz;
    angular_momentum_z = x*vy-y*vx;
    //write to the file
    writefile(t, x, y, z, vx, vy, vz,E_kinetic,E_potential,E_total,angular_momentum_x,angular_momentum_y,angular_momentum_z);
  }
  outfile.close();
  verfinish = clock();



  //velocity Verlet no EULER
  //initialize again
  t = 0.0;
  x = x0;
  y = y0;
  z = z0;
  vx = vx0;
  vy = vy0;
  vz = vz0;
  r = sqrt(x*x + y*y + z*z);
  E_kinetic = 0.;
  E_potential = 0.;
  E_total = 0.;
  angular_momentum_x = 0.;
  angular_momentum_y = 0.;
  angular_momentum_z = 0.;

  eulstart = clock();
  outfile.open("sun-earth-euler.txt");

  //time loop
  while (t <= t_fin) {
    //calculate accelerations
    ax = -Gm_sun*x/(r*r*r);
    ay = -Gm_sun*y/(r*r*r);
    az = -Gm_sun*z/(r*r*r);
    //change the values
    x += dt*vx;
    y += dt*vy;
    z += dt*vz;
    vx += dt*ax;
    vy += dt*ay;
    vz += dt*az;
    r = sqrt(x*x + y*y + z*z);
    t += dt;
    E_kinetic = 0.5*(vx*vx + vy*vy + vz*vz);
    E_potential = -Gm_sun/r;
    E_total = E_kinetic + E_potential;
    angular_momentum_x = y*vz-z*vy;
    angular_momentum_y = z*vx-x*vz;
    angular_momentum_z = x*vy-y*vx;
    //write to the file
    writefile(t, x, y, z, vx, vy, vz,E_kinetic,E_potential,E_total,angular_momentum_x,angular_momentum_y,angular_momentum_z);
  }
  outfile.close();
  eulfinish = clock();



  cout <<setprecision(4)<< "Forward Euler: " <<  double(eulfinish-eulstart)/CLOCKS_PER_SEC << "s" << endl;
  cout <<setprecision(4)<< "Velocity Verlet: " << double(verfinish-verstart)/CLOCKS_PER_SEC << "s" << endl;

  return 0;
}


void writefile(double t, double x, double y, double z, double vx, double vy, double vz, double E_kinetic, double E_potential, double E_total, double angular_momentum_x, double angular_momentum_y, double angular_momentum_z)
{
  outfile << setiosflags(ios::showpoint);
  outfile << setw(15) << setprecision(8) << t;
  outfile << setw(15) << setprecision(8) << x;
  outfile << setw(15) << setprecision(8) << y;
  outfile << setw(15) << setprecision(8) << z;
  outfile << setw(15) << setprecision(8) << vx;
  outfile << setw(15) << setprecision(8) << vy;
  outfile << setw(15) << setprecision(8) << vz;
  outfile << setw(15) << setprecision(8) << E_kinetic;
  outfile << setw(15) << setprecision(8) << E_potential;
  outfile << setw(15) << setprecision(8) << E_total;
  outfile << setw(15) << setprecision(8) << angular_momentum_x;
  outfile << setw(15) << setprecision(8) << angular_momentum_y;
  outfile << setw(15) << setprecision(8) << angular_momentum_z << endl;
}
