#include "catch.hpp"
#include "headerTest.hpp"


TEST_CASE("Testing max a(i,j)"){

    int n=3;  
    int interact = 2; // 0 for 1 electorn, 1 for two electrons, 2 buckling beam
    double w = 1; // 1 no amort 
    double RMax = 1; // RMax=1 for a-c tasks, RMax other values for following tasks
    mat a = zeros<mat>(n,n);
    vec radius(n);
    initialise(n, a, radius, RMax, w, interact);
// cout<<a;
    int p=0;
    int q=0;
    double max=0;
    //find maximum matrix element
    offdiag(a, p, q, n, max);
    REQUIRE(p==0);
    REQUIRE(q==1);
    REQUIRE(max==Approx(9));
}

TEST_CASE("Testing eigenvalues"){
    int n=4;
    double pi=acos(-1);
    int interact = 2; // 0 for 1 electron, 1 for two electrons, 2 for buckling beam
    double w = 1; // 1 no amort 
    double RMax = 1; // RMax=1 for a-c tasks, RMax=10 for unit test
    vec radius(n);
    mat A = zeros<mat>(n,n);
    mat Eigvectores = eye<mat>(n,n);
    //initialize matrix
    initialise(n, A, radius, RMax, w, interact);
    //do jacobi algorithm until convergence
    double tolerance = 1.0E-10; 
    int iterations = 0; 
    int maxiter = 2*n*n;
    int p, q;
    double maxnondiag = tolerance*2;
    while ( maxnondiag > tolerance && iterations <= maxiter)
    {
        offdiag(A, p, q, n, maxnondiag);
        Jacobi_rotate(A, Eigvectores, p, q, n);
        iterations++;
    }   
    //get eigenvalue vector
  vec NumericalValues = diagvec(A);
  sortResult (NumericalValues, Eigvectores, n);
    
    // get exact eig values/vectors
    double Step = 1.0/n;
    double DiagConst = 2.0/(Step*Step);
    double NondiagConst =  -1.0/(Step*Step);
    vec ExactValues = zeros<vec>(n);
    if(interact == 2){
        for(int i = 0; i < n; i++) {
            ExactValues(i) = DiagConst+2*NondiagConst*cos((i+1)*pi/(n+1));
        }
        cout<<ExactValues;
        REQUIRE(NumericalValues(0)==Approx(ExactValues(0)));
        REQUIRE(NumericalValues(1)==Approx(ExactValues(1)));
        REQUIRE(NumericalValues(2)==Approx(ExactValues(2)));
        REQUIRE(NumericalValues(3)==Approx(ExactValues(3)));
        REQUIRE(iterations<maxiter);
    }
    else if(interact == 0){
        REQUIRE(NumericalValues(0)==Approx(6.56863));
        REQUIRE(NumericalValues(1)==Approx(25.32055));
        REQUIRE(NumericalValues(2)==Approx(56.57082));
        REQUIRE(iterations<maxiter);
    }
}

TEST_CASE("Testing eigenvector orthogonality"){
    int n=4;
    double pi=acos(-1);
    int interact = 2; // 0 for normal, 1 for two electrons
    double w = 1; // 1 no amort 
    double RMax = 1; // RMax=1 for a-c tasks, RMax other values for following tasks
    vec radius(n);
    mat A = zeros<mat>(n,n);
    mat Eigvectores = eye<mat>(n,n);
    //initialize matrix
    initialise(n, A, radius, RMax, w, interact);
    //do jacobi algorithm until convergence
    double tolerance = 1.0E-10; 
    int iterations = 0; 
    int maxiter = 2*n*n;
    int p, q;
    double maxnondiag = tolerance*2;
    while ( maxnondiag > tolerance && iterations <= maxiter)
    {
        offdiag(A, p, q, n, maxnondiag);
        Jacobi_rotate(A, Eigvectores, p, q, n);
        iterations++;
    }   
    // cout<<Eigvectores;
    //get eigenvalue vector
  vec NumericalValues = diagvec(A);
  sortResult (NumericalValues, Eigvectores, n);
    // cout<<Eigvectores;
    // cout<<NumericalValues;
    //test eigen vector orthogonality
    //dot1=v0*v1=0
    double dot1=Eigvectores(0,0)*Eigvectores(1,0)+Eigvectores(0,1)*Eigvectores(1,1)
        +Eigvectores(0,2)*Eigvectores(1,2)+Eigvectores(0,3)*Eigvectores(1,3);
    //dot2=v0*v0=1
    double dot2=Eigvectores(0,0)*Eigvectores(0,0)+Eigvectores(0,1)*Eigvectores(0,1)
        +Eigvectores(0,2)*Eigvectores(0,2)+Eigvectores(0,3)*Eigvectores(0,3);
    REQUIRE(dot1==Approx(0.000).epsilon(0.01));
    REQUIRE(dot2==Approx(1.000));
}
