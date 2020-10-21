/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef LINEAR_EQUATION_H
#define LINEAR_EQUATION_H

#include <math.h>
#include <string>
#include "basicMatrix.h"

#define ERROR_NO_SOLUTION     2

class LinearEquation : public Utilities
{
  public:
    LinearEquation(string s="");
    LinearEquation(bool empty);
    LinearEquation(BasicMatrix a, BasicMatrix b);
    virtual ~LinearEquation() { };
    BasicMatrix matrixA, matrixB, solutionMatrix;
    unsigned int numberVariables;
    vector<double> solution;


    virtual bool Solve() {return true;};
    bool ForwardSubstitution(BasicMatrix triangular, BasicMatrix b, BasicMatrix &result);
    bool BackwardsSubstitution(BasicMatrix U, BasicMatrix Y, BasicMatrix &X);
    bool PrintSolution();
};

class Jacobi : public LinearEquation
{
  public:
    Jacobi(string s="");
    double threshold = 0.0;
    bool Solve();
};

class GaussSeidel : public LinearEquation
{
  public:
    GaussSeidel(string s="");
    double threshold = 0.0;
    bool Solve();
};


class LU : public LinearEquation
{
  public:
    LU(string s="");
    LU(bool empty) : LinearEquation(empty) {};
    LU(BasicMatrix a, BasicMatrix b) : LinearEquation (a, b) {this->Check();};
    BasicMatrix matrixLU;
    bool Decompose();
    void Check();
    void MakeL (BasicMatrix &L);
    void MakeU (BasicMatrix &U);
    bool Solve();
};

class Cholesky : public LinearEquation
{
  public:
    Cholesky(string s="");
    BasicMatrix L, LT;
    bool Decompose();
    bool Solve();
};
#endif  

