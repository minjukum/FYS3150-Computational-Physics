//tests the find_maximum function
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void find_maximum(mat, int&, int&, double&, int);

int main(int argc, char const *argv[]) {
  double Akl;
  int k,l;
  int n = 5;
  mat a = mat(n,n); a.randu();

  cout << "original matrix is: " <<endl<< a <<endl;

  find_maximum(a, k, l, Akl, n);

  cout << "The maximum: " << Akl <<endl;
  cout <<"k, l is "<< k <<", "<< l <<endl;

  return 0;
}

void find_maximum(mat A, int&k, int&l, double&Akl, int n){
  Akl = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i!=j && abs(A(i,j))>=abs(Akl)) {
        Akl = A(i,j);
        k=i;
        l=j;
      }
    }
  }
}
