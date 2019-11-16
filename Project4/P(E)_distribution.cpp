//Nov. 2019, Minju Kum
//Project 4d) code for Ising model (20x20 lattice)
//see the probability distribution of energy
//in the interval (-800, 800), bin size = 50
//at the given temperature T, after stable state reached
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

//initialization of spin matrix, E, M, weights, expectation values
void initialize(int, mat&, double&, double&, int, vec&, double, vec&);

//inline function for periodic boundary condtition
inline int periodic(int, int, int);

//Metropolis algorithm
void metropolis(int, mat&, vec, double&, double&);


int main(int argc, char const *argv[]) {

  int L = 20;
  double T;
  int randomness = 0;
  long int n = 10000000; //number of MC cycles
  long int trash_garbage = 1000000; //steady state reached after this
  double E, M;
  int hist_total = 0;

  mat spin = mat(L,L);
  vec w = vec(17);
  vec expectation = vec(5);
  vec hist = vec(32);
  hist.zeros();

  //read in temperature, # range of MC cycles
  cout << "Enter the temperature: ";
  cin >> T;
//  cout << "Would you like the initial state to be random? 1 or 0: ";
//  cin >> randomness;
//  cout << "Enter the number of Monte-Carlo cycles: ";
//  cin >> n;

  //set a filename with temperature
  ofstream histfile;
  string histname = "4d_hist_";
  string temperature = to_string(T);
  histname.append(temperature);
  histname.append(".txt");
  histfile.open(histname);

  //initialization of spin matrix, E, M, weights, expectation values
  initialize(randomness, spin, M, E, L, w, T, expectation);

  //MC cycles, n times
  for (int i = 1; i <= n; i++) {
    //metropolis algorithm
    metropolis(L, spin, w, E, M);

    //after steady state is reached
    if (i > trash_garbage) {
      //update expectation values
      expectation(0) += E;
      expectation(1) += E*E;

      //count the number of times energy appearing in the given interval
      for (int j = -16; j < 16; j++) {
        if ((j*50. <= E)&&(E < (j*50.+50.))) {
          hist(j+16) += 1;
        }
      }
      if (E == 800.) {
        hist(31) += 1;
      }
    }
  }

  //just for a check, see if (hist_total) = (n-trash_garbage)
  for (int i = 0; i < 32; i++) {
    hist_total += hist(i);
  }

  //final calculation, print the results
  //all printed values are "per one spin"
  double mean_E = expectation(0)/(n-trash_garbage);
  double var_E = expectation(1)/(n-trash_garbage) - mean_E*mean_E;
  cout << setprecision(8) << "mean E: " << mean_E/(L*L) << endl;
  cout << setprecision(8) << "variance of E: " << var_E/(L*L) << endl;

  //write the histogram data
  for (int i = -16; i < 16; i++) {
    histfile << setw(15) << setprecision(8) << i*50.;
    histfile << setw(15) << setprecision(8) << hist(i+16) << endl;
  }
  //check the total number of hist()
  cout << hist_total << " counts in histogram out of " << n - trash_garbage << " effective MC cycles" << endl;

  histfile.close();


  return 0;
}


void initialize(int randomness, mat& spin, double& M, double& E, int L, vec& w, double T, vec& expectation){
  //initialize spin matrix
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

  //set up probability(weights) array
  for (int dE = -8; dE <= 8; dE++)
    w(dE+8) = 0;
  for (int dE = -8; dE <= 8; dE+=4)
    w(dE+8) = exp(-double(dE)/T);

  //initialize expectation values
  for (int i = 0; i < 5; i++)
    expectation(i) = 0.;
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
