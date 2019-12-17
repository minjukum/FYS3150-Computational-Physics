/*
Nov.2019 Minju Kum
-This program simulates Sun-Earth system in 3d with velocity Verlet algorithm.
-Gm_sun is gravitational constant G multiplied by mass of the Sun,
 which is approximated as 4*pi*pi with the unit (1AU)^3/(1year)^2.
-Different gravitational laws were also tested with parameter beta.
 beta appears here: F = -(GMm)/(r^beta)

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

void writefile(double, double, double, double, double, double, double);

int main(int argc, char const *argv[]) {

  //variables
  long int n;
  double t_fin,dt,t;
  double r,x,y,z,vx,vy,vz;
  double Gm_sun;
  double beta;
  clock_t verstart, verfinish;

  //read in final time, number of integration points, initial conditions
  cout << "Enter the final time and number of integration points : ";
  cin >> t_fin >> n;
//  cout << "Enter the initial position x, y, z : ";
//  cin >> x >> y >> z;
//  cout << "Enter the initial velocity v_y: ";
//  cin >> vy;
  cout << "Enter beta: ";
  cin >> beta;

  x = 4.886437732350767*0.1;
  y = 8.637782835015937*0.1;
  z = -2.647417851397024*0.00001;
  vx = -1.520571286699704*365*0.01;
  vy = 8.506149010764619*365*0.001;
  vz = -1.220989097339310*365*0.0000001;

  //set constants
  dt = t_fin/n;
  Gm_sun = 4*acos(-1.0)*acos(-1.0); // = 4*pi*pi

  //velocity Verlet
  //initialize
  t = 0.0;
  r = sqrt(x*x + y*y + z*z);

  char* filename = new char[100];
  sprintf(filename, "sun-earth.txt");
  outfile.open(filename);

  verstart = clock();

  //time loop
  while (t <= t_fin) {
    //save the original positions to calculate velocities
    double x1 = x;
    double y1 = y;
    double z1 = z;
    //change the positions
    x += dt*vx - 0.5*dt*dt*Gm_sun*x/pow(r,beta+1);
    y += dt*vy - 0.5*dt*dt*Gm_sun*y/pow(r,beta+1);
    z += dt*vz - 0.5*dt*dt*Gm_sun*z/pow(r,beta+1);
    //and we also need new r
    double r_new = sqrt(x*x + y*y + z*z);
    //keep changing
    vx -= 0.5*dt*Gm_sun*(x/pow(r_new,beta+1) + x1/pow(r,beta+1));
    vy -= 0.5*dt*Gm_sun*(y/pow(r_new,beta+1) + y1/pow(r,beta+1));
    vz -= 0.5*dt*Gm_sun*(z/pow(r_new,beta+1) + z1/pow(r,beta+1));
    r = r_new;
    t += dt;
    //write to the file
    writefile(t, x, y, z, vx, vy, vz);
  }
  outfile.close();
  verfinish = clock();

  cout <<setprecision(4)<< "CPU time: " << (float(verfinish-verstart)/CLOCKS_PER_SEC) << " s" << endl;

  return 0;
}


void writefile(double t, double x, double y, double z, double vx, double vy, double vz)
{
  outfile << setiosflags(ios::showpoint);
  outfile << setw(15) << setprecision(8) << t;
  outfile << setw(15) << setprecision(8) << x;
  outfile << setw(15) << setprecision(8) << y;
  outfile << setw(15) << setprecision(8) << z;
  outfile << setw(15) << setprecision(8) << vx;
  outfile << setw(15) << setprecision(8) << vy;
  outfile << setw(15) << setprecision(8) << vz << endl;
}
