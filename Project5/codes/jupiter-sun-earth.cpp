/*
Nov.2019 Minju Kum
-This program simulates Sun-Earth-Jupiter system in 3d with velocity Verlet algorithm.
-Gm_sun is gravitational constant G multiplied by mass of the Sun,
 which is approximated as 4*pi*pi with the unit (1AU)^3/(1year)^2.
-Units of position and velocity is in AU and AU/Year.

   Sun
   ...    *     *         *         *
.///////.     *    Earth                                   Jupiter
.//^//^/.            .     *                                  @
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
ofstream jupiterfile;

void writefile(double, double, double, double, double, double, double);
void writefile2(double, double, double, double, double, double, double);

int main(int argc, char const *argv[]) {

  //variables
  long int n; //number of integration points
  double t_fin,dt,t;
  double x,y,z,vx,vy,vz;
  double xj,yj,zj,vxj,vyj,vzj;
  double r,rj,rej;
  double Gm_sun,Gm_jup,Gm_ear;
  double beta = 2.;
  double scaler;
  clock_t verstart, verfinish;

  //read in final time, number of integration points, initial conditions
  cout << "Enter the final time and number of integration points : ";
  cin >> t_fin >> n;
  cout << "How heavy is Jupiter? : ";
  cin >> scaler;

  //initial ctn of the Earth
  x = 4.886437732350767*0.1;
  y = 8.637782835015937*0.1;
  z = -2.647417851397024*0.00001;
  vx = -1.520571286699704*365*0.01;
  vy = 8.506149010764619*365*0.001;
  vz = -1.220989097339310*365*0.0000001;

  //initial ctn of the Jupiter
  xj = 2.324521172578399*0.1;
  yj = -5.228871634333293;
  zj = 1.648517448525890*0.01;
  vxj = 7.446401219049974*365*0.001;
  vyj = 6.947253100564677*365*0.0001;
  vzj = -1.694754869444451*365*0.0001;

  //set constants
  dt = t_fin/n;
  Gm_sun = 4*acos(-1.0)*acos(-1.0); // = 4*pi*pi
  Gm_jup = Gm_sun*954.79194*pow(10,-6)*scaler;
  Gm_ear = Gm_sun*3.00348959632*pow(10,-6);

  //velocity Verlet
  //initialize
  t = 0.0;
  r = sqrt(x*x + y*y + z*z);
  rj = sqrt(xj*xj + yj*yj + zj*zj);
  rej = sqrt((x-xj)*(x-xj)+(y-yj)*(y-yj)+(z-zj)*(z-zj));


  char* filename = new char[100];
  sprintf(filename,"jupiter-sun-earth_earth.%f.txt",scaler);
  outfile.open(filename);
  jupiterfile.open("jupiter-sun-earth_jupiter.txt");

  verstart = clock();

  //time loop
  while (t <= t_fin) {
    //save the original positions
    double x1 = x;
    double y1 = y;
    double z1 = z;
    double xj1 = xj;
    double yj1 = yj;
    double zj1 = zj;

    //change the positions
    x += dt*vx - 0.5*dt*dt*(Gm_sun*x/pow(r,beta+1) + Gm_jup*(x-xj)/pow(rej,beta+1));
    y += dt*vy - 0.5*dt*dt*(Gm_sun*y/pow(r,beta+1) + Gm_jup*(y-yj)/pow(rej,beta+1));
    z += dt*vz - 0.5*dt*dt*(Gm_sun*z/pow(r,beta+1) + Gm_jup*(z-zj)/pow(rej,beta+1));
    xj += dt*vxj - 0.5*dt*dt*(Gm_sun*xj/pow(rj,beta+1) + Gm_ear*(xj-x1)/pow(rej,beta+1));
    yj += dt*vyj - 0.5*dt*dt*(Gm_sun*yj/pow(rj,beta+1) + Gm_ear*(yj-y1)/pow(rej,beta+1));
    zj += dt*vzj - 0.5*dt*dt*(Gm_sun*zj/pow(rj,beta+1) + Gm_ear*(zj-z1)/pow(rej,beta+1));

    //and we also need new r
    double r_new = sqrt(x*x + y*y + z*z);
    double rj_new = sqrt(xj*xj + yj*yj + zj*zj);
    double rej_new = sqrt((x-xj)*(x-xj)+(y-yj)*(y-yj)+(z-zj)*(z-zj));

    //keep changing
    vx -= 0.5*dt * (Gm_sun*(x/pow(r_new,beta+1) + x1/pow(r,beta+1)) + Gm_jup*((x-xj)/pow(rej_new,beta+1) + (x1-xj1)/pow(rej,beta+1)));
    vy -= 0.5*dt * (Gm_sun*(y/pow(r_new,beta+1) + y1/pow(r,beta+1)) + Gm_jup*((y-yj)/pow(rej_new,beta+1) + (y1-yj1)/pow(rej,beta+1)));
    vz -= 0.5*dt * (Gm_sun*(z/pow(r_new,beta+1) + z1/pow(r,beta+1)) + Gm_jup*((z-zj)/pow(rej_new,beta+1) + (z1-zj1)/pow(rej,beta+1)));
    vxj -= 0.5*dt * (Gm_sun*(xj/pow(rj_new,beta+1) + xj1/pow(rj,beta+1)) + Gm_ear*((xj-x)/pow(rej_new,beta+1) + (xj1-x1)/pow(rej,beta+1)));
    vyj -= 0.5*dt * (Gm_sun*(yj/pow(rj_new,beta+1) + yj1/pow(rj,beta+1)) + Gm_ear*((yj-y)/pow(rej_new,beta+1) + (yj1-y1)/pow(rej,beta+1)));
    vzj -= 0.5*dt * (Gm_sun*(zj/pow(rj_new,beta+1) + zj1/pow(rj,beta+1)) + Gm_ear*((zj-z)/pow(rej_new,beta+1) + (zj1-z1)/pow(rej,beta+1)));

    r = r_new;
    rj = rj_new;
    rej = rej_new;

    t += dt;

    //write to the file
    writefile(t, x, y, z, vx, vy, vz);
    writefile2(t, xj, yj, zj, vxj, vyj, vzj);
  }
  outfile.close();
  jupiterfile.close();
  verfinish = clock();

  cout <<setprecision(4)<< "CPU time: " << double(verfinish-verstart)/CLOCKS_PER_SEC << "s" << endl;

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

void writefile2(double t, double x, double y, double z, double vx, double vy, double vz)
{
  jupiterfile << setiosflags(ios::showpoint);
  jupiterfile << setw(15) << setprecision(8) << t;
  jupiterfile << setw(15) << setprecision(8) << x;
  jupiterfile << setw(15) << setprecision(8) << y;
  jupiterfile << setw(15) << setprecision(8) << z;
  jupiterfile << setw(15) << setprecision(8) << vx;
  jupiterfile << setw(15) << setprecision(8) << vy;
  jupiterfile << setw(15) << setprecision(8) << vz << endl;
}
