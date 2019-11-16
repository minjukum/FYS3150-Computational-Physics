//Nov. 2019, Minju Kum
//Project 4c) code for Ising model (20x20 lattice)
//computes the mean energy E, mean magnetization |M|
//and see how many MC cycles are needed to reach the equilibrium
//at the given temperature T
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
void metropolis(int, mat&, vec, double&, double&, long int&);


int main(int argc, char const *argv[]) {

  int L = 20;
//  double T_init, T_fin, dT;
  double T;
  int randomness;
  int n_init, n_fin;
  long int n;
  long int n_accepted;
  double E, M;


  mat spin = mat(L,L);
  vec w = vec(17);
  vec expectation = vec(5);

  //read in temperature, randomness, # range of MC cycles
  cout << "Enter the temperature: ";
  cin >> T;
  cout << "Would you like the initial state to be random? 1 or 0: ";
  cin >> randomness;
  cout << "Enter the number range of Monte-Carlo cycles in log scale: ";
  cin >> n_init >> n_fin;

  //set a filename with temperature
  ofstream outfile;
  string filename = "4c_";
  string temperature = to_string(T);
  filename.append(temperature);
  if (randomness == 1) {
    string ran = "ran";
    filename.append(ran);
  }
  filename.append(".txt");
  outfile.open(filename);

  //n loop---------------------------------------
  //span over the given number range of MC cycles
  for (int exponent = n_init; exponent <= n_fin; exponent++) {
    //set the number of MC cycles, initialize n_accpeted
    n = pow(10,exponent);
    n_accepted = 0;

    //initialize spin matrix, magnetization, energy
    if (randomness == 1) {
      spin.randu();
      spin.clean(0.5);
      spin.for_each([](mat::elem_type& val){val -= 1;});
      spin.clean(0.5);
      spin.replace(0.0,1.0); //random state
    }
    else spin.ones(); //ground state

    //initial calculation of magnetization and energy
    M = E = 0.;
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
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
    for (int dE = -8; dE <= 8; dE+=4)
      w(dE+8) = exp(-double(dE)/T);

    //initialize expectation values
    for (int i = 0; i < 5; i++)
      expectation(i) = 0.;

    //MC cycles, n times
    for (int i = 0; i < n; i++) {
      //metropolis algorithm
      metropolis(L, spin, w, E, M, n_accepted);
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
    outfile << setw(15) << n;
    outfile << setw(15) << setprecision(8) << mean_E/(L*L); // mean E
    outfile << setw(15) << setprecision(8) << var_E/(T*T*L*L);  // C_v
    outfile << setw(15) << setprecision(8) << mean_M/(L*L); // mean M
    outfile << setw(15) << setprecision(8) << var_M/(T*L*L); // susceptibility
    outfile << setw(15) << setprecision(8) << mean_fabsM/(L*L); // mean |M|
    outfile << setw(15) << setprecision(8) << double(n_accepted)/(L*L);//# of accepted configurations
    outfile << setw(15) << setprecision(8) << sqrt(var_E/n)/(L*L);  // variance of mean_E
    outfile << setw(15) << setprecision(8) << sqrt(var_M/n)/(L*L) <<endl;  // variance of mean_fabsM

  } //end of n loop-------------------

  outfile.close();


  return 0;
}


void metropolis(int L, mat& spin, vec w, double& E, double& M, long int& n_accepted){
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
        n_accepted += 1;
      }
    }
  }
}

inline int periodic(int i, int limit, int add){
  return (i+limit+add)%limit;
}
