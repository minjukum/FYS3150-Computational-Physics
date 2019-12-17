#ifndef MYPLANET_H
#define MYPLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class planet
{
public:
    // Properties
    double mass;
    double position[3];
    double velocity[3];
    double position_old[3];
    double potential;
    double kinetic;

    // Initializers
    planet();
    planet(double M,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    void save_position();
    double distance(planet otherPlanet);
    double distance_old(planet otherPlanet);

};

#endif // MYPLANET_H
