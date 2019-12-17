#ifndef MYSOLVER_H
#define MYSOLVER_H
#define _USE_MATH_DEFINES
#include "myplanet.h"
#include <vector>
#include <fstream>
using std::vector;

class solver
{
public:
    friend class planet;

    // properties
    double radius,total_mass,G;
    int total_planets;
    vector<planet> all_planets;
    double totalKinetic;
    double totalPotential;

    // constants

    // initializers
    solver();

    // functions
    void add(planet newplanet);
    void writefile(std::ofstream &outfile, double time, int number);
    void VelocityVerlet(double t_fin,long int integration_points,int number,int is_sun_moving);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);

};

#endif // MYSOLVER_H
