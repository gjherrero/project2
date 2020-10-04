 // 
// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <armadillo>
// #include "catch.hpp"
// #define CATCH_CONFIG_RUNNER
// using namespace  std;
// using namespace  arma;
// void initialise(int n,mat &A);
// void offdiag(mat A, int &p, int &q, int n,double &max);
// void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n);
// void sortResult ( vec &a, mat &R,int n);


#include "headerTest.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
double pi= acos(-1);
using namespace  std;
using namespace  arma;
// // Begin of main program   


int main(int argc, char* argv[])
{
  // loop over different max R and Mat sizes. Commented for interact==2"
  // string filename;
  // ofstream myfile;
  // for (int RMax = 10; RMax < 1000; RMax+=100)  {
  //   filename = "RMax";
  //   filename.append(to_string(RMax));
  //   filename.append(".txt");
  //   myfile.open(filename);
  int factor = 1;
  //   for(int Dim=factor*10; Dim < factor*110; Dim+=factor*10){
      int Dim = 100;
      int RMax = factor*10;
      int interact = 0; // 0 for 1 electron, 1 for two electrons //2 for Bucling beam
      double w = 1; // 1 no amort 
      mat Hamiltonian = zeros<mat>(Dim,Dim);
      mat Eigvectores = eye<mat>(Dim,Dim);
      vec radius(Dim);
      initialise(Dim, Hamiltonian, radius, RMax, w, interact);
      // cout<<Hamiltonian;
      // armadillo values and vectors
      vec eigval;
      mat eigvec;
      mat Hamiltemp = Hamiltonian;
        clock_t t0, t1, t2;
        t0=clock();
        eig_sym( eigval, eigvec, Hamiltemp);
        t1=clock();

      // own code for eig values/vectors only used for interact==2
      if (interact ==2){      
        double tolerance = 1.0E-10; 
        int iterations = 0; 
        int maxiter = 2*Dim*Dim;
        int p, q;
        double maxnondiag = tolerance*2;
        while ( maxnondiag > tolerance && iterations <= maxiter)
        {
          offdiag(Hamiltonian, p, q, Dim, maxnondiag);
          Jacobi_rotate(Hamiltonian, Eigvectores, p, q, Dim);
          iterations++;
        }   
        // sort results and compare
        vec NumericalValues = diagvec(Hamiltonian);
        sortResult (NumericalValues, Eigvectores, Dim);
        t2=clock();
            // get exact eig values/vectors
        mat ExactVectors = zeros<mat>(Dim,Dim);
        vec ExactValues = zeros<vec>(Dim);
        double h = RMax/ double(Dim);
        for(int i = 0; i < Dim; i++) {
          ExactValues(i) = (2/(h*h))+2*(-1/(h*h))*cos((i+1)*pi/(Dim+1));
          for(int j = 0; j < Dim; j++) {
            ExactVectors(i,j) = sin((i+1)*(j+1)*pi/(Dim+1));
          }
        }        
        cout<<Dim<<"  "<<iterations<<"  "<<((double)t2-(double)t1)/CLOCKS_PER_SEC<<"  "<<((double)t1-(double)t0)/CLOCKS_PER_SEC<<";"<<endl;
      }
        
          cout<<eigval<<endl;
          cout<<eigvec<<endl; 
    
      return 0;
      // printing results to file for interact==0 and interact ==1 comment if interact==2
  //     for(int i=0;i<4;i++){
  //         // myfile << setw(15) << setprecision(8) << NumericalValues[i];
  //         myfile << setw(15) << setprecision(8) << eigval[i];
  //     }
  //     myfile << endl;
  //     cout<< Dim<<endl;
  //  }
  //   myfile.close();
  // }
}