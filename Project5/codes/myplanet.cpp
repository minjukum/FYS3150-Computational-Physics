#include "myplanet.h"

planet::planet()
{
    mass = 1.;
    position[0] = 1.;
    position[1] = 0.;
    position[2] = 0.;
    position_old[0] = 0.;
    position_old[1] = 0.;
    position_old[2] = 0.;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    potential = 0.;
    kinetic = 0.;
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    position_old[0] = 0.;
    position_old[1] = 0.;
    position_old[2] = 0.;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
}

void planet::save_position()
{
  for (int i = 0; i < 3; i++) {
    this->position_old[i] = this->position[i];
  }
}

double planet::distance(planet otherPlanet)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = otherPlanet.position[0];
    y2 = otherPlanet.position[1];
    z2 = otherPlanet.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
 }

double planet::distance_old(planet otherPlanet)
{
  double x1,y1,z1,x2,y2,z2,xx,yy,zz;

  x1 = this->position_old[0];
  y1 = this->position_old[1];
  z1 = this->position_old[2];

  x2 = otherPlanet.position_old[0];
  y2 = otherPlanet.position_old[1];
  z2 = otherPlanet.position_old[2];

  xx = x1-x2;
  yy = y1-y2;
  zz = z1-z2;

  return sqrt(xx*xx + yy*yy + zz*zz);
}
