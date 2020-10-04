#include "headerTest.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "catch.hpp"
using namespace  std;
using namespace  arma;

void initialise(int n, mat &A,vec &r, double RMax, double w, int interact){ //Add Rmax and also add vector (as jacobi from teacher) initialise with included V_r as a funcion of vector positions      
  double h = RMax/ double(n);    
  //initialize r values
  r(0)=h;
  for (int i=1; i<n ;i++){
    r(i)=r(i-1)+h;
  }

  //initialize matrix and vector
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      if(i==j && interact==0){
        A(i,j)=2/(h*h)+w*w*r(i)*r(i);
      }
      else if (i==j && interact==1){
        A(i,j)=2/(h*h)+w*w*r(i)*r(i)+1/r(i);
      }
      else if (i==j && interact==2){
        A(i,j)=2/(h*h);
      }
      else if (i==j+1 or i==j-1){
        A(i,j)=-1/(h*h);
      } 
      else{
        A(i,j)=0;
      }         
    }
  }
}
void sortResult( vec &a, mat &R,int n){
  uvec indices = sort_index(a);
  mat R_temp = R;
  for ( int i = 0; i < n; i++){
    for ( int j = 0; j < n; j++){
      R(i,j)=R_temp(i,indices(j));
    }
  }
  a = sort(a);
}

//  the offdiag function, using Armadillo
void offdiag(mat A, int &p, int &q, int n,double &max){
   max = 0.0;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j) //starts i+1, only searchs upper right side of matrix
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           { 
              max = aij;  p = i; q = j;
           }
       }
   }
}

void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n )
{
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
  return;
} // end of function jacobi_rotate