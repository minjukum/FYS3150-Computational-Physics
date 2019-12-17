/*
Nov.2019 Minju Kum
-This program simulates Sun-Earth-Jupiter system in 3d with velocity Verlet algorithm
-Object oriented
-All mass is relative to the Sun(/M_sun)
-Units of position and velocity is in AU and AU/Year

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
#include "myplanet.h"
#include "mysolver.h"

using namespace std;

int main(int argc, char const *argv[]) {

  //variables
  long int integration_points;
  int is_sun_moving;
  int number = 3; //print results for all three bodies
  double t_fin;

  //read in final time, number of integration points
  cout << "Enter the final time and number of integration points: ";
  cin >> t_fin >> integration_points;
  //check if the Sun is also moving
  cout << "Is the Sun moving? (1/0): ";
  cin >> is_sun_moving;

  //create planets
  planet Sun(1,0,0,0,0,0,0);
  planet Earth(3.00348959632*pow(10,-6),1.,0.,0.,0.,2*acos(-1),0.);
  planet Jupiter(954.79194*pow(10,-6)*1000,5.203,0.,0.,0.,2.76,0.);

  //create system and add planets
  solver SolarSystem;
  SolarSystem.add(Sun);
  SolarSystem.add(Earth);
  SolarSystem.add(Jupiter);

  //velocity verlet
  SolarSystem.VelocityVerlet(t_fin,integration_points,number,is_sun_moving);

  return 0;
}
