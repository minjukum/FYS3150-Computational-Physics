//Sep.2019 written by Minju Kum
//notes:
//1.these functions are written using armadillo library,
//  which stores matrices in column major ordering.
//2.some for loops that inlcude abs(i-j) are defined with int i,j.
//  others are defined with uword i,j

#include "jacobi.h"

//initialize vector and matrix, buckling beam
void initialize_beam(mat& A, vec& rho, int n, double h){

  //initialize vector
  rho(0) = 0;
  rho(n+1) = 1;
  for (int i = 1; i < n+1; i++) {
    rho(i) = rho(0) + i*h;
  }

  //initialize matrix
  double h_squared = h*h;
  for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        if (i == j)
          A(i,j) = 2/h_squared;
        else if (abs(i-j) == 1)
          A(i,j) = -1/h_squared;
        else
          A(i,j) = 0.0;
      }
  }
}

//initialize vector and matrix, harmonic potential, 1 electron
void initialize_harmonic_1el(mat& A, vec& rho, int n, double h, double rho_end){

  //initialize vector
  rho(0) = 0;
  rho(n+1) = rho_end;
  for (int i = 1; i < n+1; i++) {
    rho(i) = rho(0) + i*h;
  }

  //initialize matrix
  double h_squared = h*h;
  for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        if (i == j)
          A(i,j) = 2/h_squared + rho(i+1)*rho(i+1);
        else if (abs(i-j) == 1)
          A(i,j) = -1/h_squared;
        else
          A(i,j) = 0.0;
      }
  }
}

//initialize vector and matrix, harmonic potential, 2 electrons
void initialize_harmonic_2el(mat& A, vec& rho, int n, double h, double w, double interaction, double rho_end){

  //initialize vector
  rho(0) = 0;
  rho(n+1) = rho_end;
  for (uword i = 1; i < n+1; i++) {
    rho(i) = rho(0) + i*h;
  }

  //initialize matrix, interaction is either 1 or 0
  double h_squared = h*h;
  double w_squared = w*w;
  for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        if (i == j)
          A(i,j) = 2/h_squared + w_squared*rho(i+1)*rho(i+1) + interaction/rho(i+1);
        else if (abs(i-j) == 1)
          A(i,j) = -1/h_squared;
        else
          A(i,j) = 0.0;
      }
  }

}

//find the greatest nondiagonal
void find_maximum(mat A, int&k, int&l, double&Akl, int n){
  Akl = 0;
  for (uword j = 0; j < n; j++) {
    for (uword i = 0; i < n; i++) {
      if (i!=j && abs(A(i,j))>=abs(Akl)) {
        Akl = A(i,j);
        k=i;
        l=j;
      }
    }
  }
}

//jacobi method, diagonalize A and compute eigenvectors
void jacobi(mat& A, mat&V, int n, double Akl, int k, int l){

  double tau, tan, cos, sin;
  double Akk, All, Aik, Ail;
  double Vik, Vil;

  //assign trigonometric ftn values
  Akk = A(k,k);
  All = A(l,l);
  tau = (All-Akk)/(2*Akl);
  if (tau >= 0)
    tan = -tau - sqrt(1+tau*tau);
  else
    tan = -tau + sqrt(1+tau*tau);
  cos = 1 / sqrt(1+tan*tan);
  sin = cos * tan;

  //change the matrix elements of A
  A(k,k) = Akk*cos*cos - 2*Akl*cos*sin + All*sin*sin;
  A(l,l) = All*cos*cos + 2*Akl*cos*sin + Akk*sin*sin;
  A(k,l) = 0;
  A(l,k) = 0;
  for (uword i = 0; i < n; i++) {
    if (i!=k && i!=l) {
      Aik = A(i,k);
      Ail = A(i,l);
      A(i,k) = Aik*cos - Ail*sin;
      A(k,i) = A(i,k);
      A(i,l) = Ail*cos + Aik*sin;
      A(l,i) = A(i,l);
    }
  }

  //change the matrix elements of V
  for (uword i = 0; i < n; i++) {
    Vik = V(i,k);
    Vil = V(i,l);
    V(i,k) = cos*Vik - sin*Vil;
    V(i,l) = cos*Vil + sin*Vik;
  }

}

//if n>=10, print the first 3 eigenvalues and save 3 eigenvectors
//else print all eigenvalues and save all eigenvectors
void printout(int n, mat A, mat V){
  if (n>=10) {
    //3 eigenvalues
    vec d = A.diag();
    vec d_sorted = sort(d);
    vec eigvals3 = d_sorted.head(3);

    //3 eigenvectors
    mat eigvecs3 = mat(n,3);
    for (uword j = 0; j < 3; j++) {
      for (uword i = 0; i < n; i++) {
        if (d(i) == d_sorted(j)) {
          eigvecs3.col(j) = V.col(i);
        }
      }
    }
    //print, save
    cout << "The first 3 eigenvalues are: " <<endl<< eigvals3 <<endl;
    eigvecs3.save("eigvecs3.txt", raw_ascii);
  }

  else {
    cout << "Diagonalized A is: " <<endl<< A << endl;
    V.save("eigvecs3.txt", raw_ascii);
  }

}
