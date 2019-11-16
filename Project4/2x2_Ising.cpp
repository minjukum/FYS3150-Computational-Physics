//Nov. 2019, Minju Kum
//Project 4b) code for Ising model (2x2 lattice)
//computes the mean energy E, mean magnetization |M|,
//specific heat C_v, susceptibility as functions of T
//use periodic boundary condition, set k,J = 1
//***compile with C++11***

#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

//inline function for periodic boundary condtition
inline int periodic(int, int, int);
//Metropolis algorithm
void metropolis(int, mat&, vec, double&, double&);


int main(int argc, char const *argv[]) {

  int L = 2;
//  double T_init, T_fin, dT;
  double T;
  int n_init, n_fin;
  long int n;
  double E, M;

  mat spin = mat(L,L);
  vec w = vec(17);
  vec expectation = vec(5);

  ofstream outfile;

  outfile.open("4b.txt");

  //read in temperature, # range of MC cycles
  cout << "Enter the temperature: ";
  cin >> T;
  cout << "Enter the number range of Monte-Carlo cycles in log scale: ";
  cin >> n_init >> n_fin;

  //n loop---------------------------------------
  //span over the given number range of MC cycles
  for (int exponent = n_init; exponent <= n_fin; exponent++) {
    //set the number of MC cycles
    n = pow(10,exponent);

    //initialize spin matrix, magnetization, energy
    E = M = 0.;
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        spin(j,i) = 1.0; //ground state
        M += double(spin(j,i));
      }
    }
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        E -= double(spin(j,i)*(spin(periodic(j,L,-1),i)+spin(j,periodic(i,L,-1))));
      }
    }

    //set up probability array
    for (int dE = -8; dE <= 8; dE++)
      w(dE+8) = 0;
    for (int dE = -8; dE <= 8; dE+=8)
      w(dE+8) = exp(-double(dE)/T);

    //initialize array for expectation values
    for (int i = 0; i < 5; i++)
      expectation(i) = 0.;

    //MC cycles, n times
    for (int i = 0; i < n; i++) {
      //metropolis algorithm
      metropolis(L, spin, w, E, M);
      //update expectation values
      expectation(0) += E;
      expectation(1) += E*E;
      expectation(2) += M;
      expectation(3) += M*M;
      expectation(4) += fabs(M);
    }

    //final calculation, write the results to the file
    //all printed values are "per one spin"
    double mean_E = expectation(0)/n;
    double var_E = expectation(1)/n - mean_E*mean_E;
    double mean_M = expectation(2)/n;
    double mean_fabsM = expectation(4)/n;
    double var_M = expectation(3)/n - mean_fabsM*mean_fabsM;
    outfile << setiosflags(ios::showpoint);
    outfile << setw(15) << setprecision(8) << n;
    outfile << setw(15) << setprecision(8) << mean_E/(L*L); // mean E
    outfile << setw(15) << setprecision(8) << var_E/(T*T*L*L);  // C_v
    outfile << setw(15) << setprecision(8) << mean_M/(L*L); // mean M
    outfile << setw(15) << setprecision(8) << var_M/(T*L*L); // susceptibility
    outfile << setw(15) << setprecision(8) << mean_fabsM/(L*L) << endl; // mean |M|

  } //end of n loop-------------------

  outfile.close();

  //analytical calculations, T = 1.0
  double z = 2*exp(8.)+12+2*exp(-8.);
  double an_mean_E = (-16*exp(8.)+16*exp(-8.))/z;
  double an_mean_absM = (8*exp(8.)+16)/z;
  double an_C_v = (64*exp(8.)+128*exp(-8.)+64*exp(8.))/z - an_mean_E*an_mean_E;
  double an_susceptibility = (32*exp(8.)+32)/z - an_mean_absM*an_mean_absM;

  cout << "analytical values are: " << endl;
  cout << setw(15) << setprecision(8) << an_mean_E/(L*L);
  cout << setw(15) << setprecision(8) << an_C_v/(L*L);
  cout << setw(15) << setprecision(8) << 0.;
  cout << setw(15) << setprecision(8) << an_susceptibility/(L*L);
  cout << setw(15) << setprecision(8) << an_mean_absM/(L*L) << endl;

  return 0;
}


void metropolis(int L, mat& spin, vec w, double& E, double& M){
  //set up random number generator, seed given as number of ticks
  //uniform distribution in the range of (0,1)
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  uniform_real_distribution<double> distribution(0.,1.);

  //repeat for # of spins
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      //pick one random spin
      int x = int(distribution(generator)*double(L));
      int y = int(distribution(generator)*double(L));

      //calculate energy difference
      int dE = 2*spin(y,x)*(spin(y,periodic(x,L,-1))+spin(periodic(y,L,-1),x)+spin(y,periodic(x,L,1))+spin(periodic(y,L,1),x));

      //Metropolis test
      if (distribution(generator) <= w(dE+8)) {
        //flip the spin
        spin(y,x) *= -1;
        //update the energy and magnetization
        E += double(dE);
        M += double(2*spin(y,x));
      }
    }
  }
}

inline int periodic(int i, int limit, int add){
  return (i+limit+add)%limit;
}
