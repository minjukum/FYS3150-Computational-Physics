#include "mysolver.h"
#include "myplanet.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include "time.h"

using namespace std;

solver::solver()
{
    total_planets = 0;
    radius = 100;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
}

void solver::add(planet newplanet)
{
    total_planets += 1;
    total_mass += newplanet.mass;
    all_planets.push_back(newplanet);
}

void solver::writefile(std::ofstream &outfile, double time, int number)
{   // Writes position and velocity of all planets to one file
  for (int i = 0; i < number; i++) {
    planet &current = all_planets[i];
    outfile << setiosflags(ios::showpoint);
    outfile << setw(15) << setprecision(8) << time;
    outfile << setw(15) << setprecision(8) << i;
    for (int i = 0; i < 3; i++) outfile << setw(15) << setprecision(8) << current.position[i];
    for (int i = 0; i < 3; i++) outfile << setw(15) << setprecision(8) << current.velocity[i];
    outfile << std::endl;
  }
}

void solver::VelocityVerlet(double t_fin,long int integration_points,int number,int is_sun_moving)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    // Define dt, initialize time
    double dt = t_fin/double(integration_points);
    double time = 0.0;

    // Set up arrays (all initialized with 0.0)
    double **acceleration_old = setup_matrix(total_planets,3);
    double **acceleration_new = setup_matrix(total_planets,3);

    // Create a file
    char *filename = new char[100];
    sprintf(filename,"%d_planets_%li.5f.txt",total_planets,integration_points);
    std::ofstream outfile(filename);

    // Write initial values to file
    writefile(outfile,time,number);

    // Set up clock to measure the time usage
    clock_t start,finish;
    start = clock();

    // PLANET CALCULATIONS
    // Loop over time
    time += dt;
    while(time <= t_fin){

      // initialize accelerations
      for(int i = 0; i < total_planets; i++){
          for(int j = 0; j < 3; j++){
              acceleration_old[i][j] = 0.0;
              acceleration_new[i][j] = 0.0;
          }
      }

      // 1. save the old positions
      for (int i = 0; i < total_planets; i++) {
        planet &current = all_planets[i];
        current.save_position();
      }

      //if sun is moving, start from the sun
      //else start from the first "planet"
      int startpoint;
      if (is_sun_moving == 1) startpoint = 0;
      else startpoint = 1;

      // 2. change the positions
      for (int p1 = startpoint; p1 < total_planets; p1++) {
        //current planet
        planet& current = all_planets[p1];
        //for each direction
        for (int i = 0; i < 3; i++) {
          //calculate acceleration from all the other planets
          for (int p2 = 0; p2 < total_planets; p2++) {
            //planets should be different
            if (p1 != p2) {
              planet& other = all_planets[p2];
              acceleration_old[p1][i] += G*other.mass*(other.position_old[i]-current.position_old[i])/pow(current.distance_old(other),3);
            }
          }
          //change position
          current.position[i] += dt*current.velocity[i] + 0.5*dt*dt*acceleration_old[p1][i];
        }
      }

      // 3. change the velocities
      for (int p1 = startpoint; p1 < total_planets; p1++) {
        //current planet
        planet& current = all_planets[p1];
        //for each direction
        for (int i = 0; i < 3; i++) {
          //calculate new acceleration from all the other planets
          for (int p2 = 0; p2 < total_planets; p2++) {
            //planets should be different
            if (p1 != p2) {
              planet& other = all_planets[p2];
              acceleration_new[p1][i] += G*other.mass*(other.position[i]-current.position[i])/pow(current.distance(other),3);
            }
          }
          //change velocity
          current.velocity[i] += 0.5*dt*(acceleration_new[p1][i]+acceleration_old[p1][i]);
        }
      }

      //write time,planet#,position,velocity to the file
      writefile(outfile, time, number);

      time += dt;
    }//end of time loop

    // Stop clock and print out time usage
    finish = clock();
    std::cout << "Total time = " << "\t" << (float(finish - start)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    std::cout << "One time step = " << "\t" << (float(finish - start)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

    // Close files
    outfile.close();

    // Clear memory
    delete_matrix(acceleration_old);
    delete_matrix(acceleration_new);
}

double ** solver::setup_matrix(int height,int width)
{   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void solver::delete_matrix(double **matrix)
{   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}
