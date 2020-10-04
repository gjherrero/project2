#ifndef headerTest_H
#define	headerTest_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <string>
#include "catch.hpp"
using namespace  std;
using namespace  arma;
void initialise(int n, mat &A,vec &r, double RMax, double w, int interact);
void offdiag(mat A, int &p, int &q, int n,double &max);
void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n);
void sortResult ( vec &a, mat &R,int n);

#endif /* JACOBI_H */