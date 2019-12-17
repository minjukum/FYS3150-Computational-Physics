/*
Nov.2019 Minju Kum
-This program simulates Solar system(Sun + 9 planets) in 3d with velocity Verlet algorithm
-Object oriented
-All mass is relative to the Sun(/M_sun)
-Units of position and velocity is in AU and AU/Year

   Sun              *                                 *
   ...    *               *         *
.///////.   M V  E  M           Jupiter     Saturn     Uranus  Neptune     Pluto
.//^//^/.   . .  .  . *            @         -O-         ø        o          .
  .///.   *       *                         *
   *                         *
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
  int is_sun_moving = 1; // Sun is NOT fixed at the center
  int number = 10; //print results for the sun and all 9 planets
  double t_fin;

  //read in final time, number of integration points
  cout << "Enter the final time and number of integration points: ";
  cin >> t_fin >> integration_points;

  //create planets
  planet Sun(1,-3.471113008766750*0.001,7.509468464928616*0.001,1.400821456788236*0.00001,-8.436791145479849*0.000001,-1.570025842052659*0.000001,2.315735577565323*0.0000001);
  planet Mercury(0.16601*pow(10,-6),-1.580774834158366*0.1,2.867526968430482*0.1,3.701504029255228*0.01,-3.028542563909174*365*0.01,-1.255664049886210*365*0.01,1.751771708721328*365*0.001);
  planet Venus(2.4478383*pow(10,-6),3.856012308211278*0.1,-6.071897313025226*0.1,-3.087275907539941*0.01,1.694691857664187*365*0.01,1.074641681932823*365*0.01,-8.307368202106533*365*0.0001);
  planet Earth(3.00348959632*pow(10,-6),4.886437732350767*0.1,8.637782835015937*0.1,-2.647417851397024*0.00001,-1.520571286699704*365*0.01,8.506149010764619*365*0.001,-1.220989097339310*365*0.0000001);
  planet Mars(0.3227151*pow(10,-6),-1.570845751722628,-4.310957204660689*0.1,2.927988829706918*0.01,4.284396107751554*365*0.001,-1.228100710648438*365*0.01,-3.624088636996528*365*0.0001);
  planet Jupiter(954.79194*pow(10,-6),2.324521172578399*0.1,-5.228871634333293,1.648517448525890*0.01,7.446401219049974*365*0.001,6.947253100564677*365*0.0001,-1.694754869444451*365*0.0001);
  planet Saturn(285.8860*pow(10,-6),3.603260982461424,-9.360263674129357,1.931247330031657*0.01,4.897563494877299*365*0.001,1.987750328038748*365*0.001,-2.295538280437610*365*0.0001);
  planet Uranus(43.66244*pow(10,-6),1.631048148335351*10,1.126774695829079*10,-1.694556380875939*0.1,-2.264414045933132*365*0.001,3.052705795597817*365*0.001,4.067378007270122*365*0.00001);
  planet Neptune(51.51389*pow(10,-6),2.921355770038784*10,-6.480050496913639,-5.398112840772882*0.1,6.589508034205666*365*0.0001,3.083299016976726*365*0.001,-7.868028827321828*365*0.00001);
  planet Pluto(0.00655*pow(10,-6),1.285693695608717*10,-3.138327842809804*10,-3.607882172011476*0.1,2.975866812996480*365*0.001,5.224578930599930*365*0.0001,-9.167097338859692*365*0.0001);

  //create system and add planets
  solver SolarSystem;
  SolarSystem.add(Sun);
  SolarSystem.add(Mercury);
  SolarSystem.add(Venus);
  SolarSystem.add(Earth);
  SolarSystem.add(Mars);
  SolarSystem.add(Jupiter);
  SolarSystem.add(Saturn);
  SolarSystem.add(Uranus);
  SolarSystem.add(Neptune);
  SolarSystem.add(Pluto);

  //velocity verlet
  SolarSystem.VelocityVerlet(t_fin,integration_points,number,is_sun_moving);

  return 0;
}
