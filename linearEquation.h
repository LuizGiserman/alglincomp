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
    LinearEquation();
    virtual ~LinearEquation() { };
    BasicMatrix matrixA, matrixB;
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
    Jacobi();
    double threshold = 0.0;
    bool Solve();
    double Residue(BasicMatrix auxiliar, BasicMatrix x0);
};

class GaussSeidel : public LinearEquation
{
  public:
    GaussSeidel();
    bool Solve();
};


class LU : public LinearEquation
{
  public:
    LU();
    BasicMatrix matrixLU;
    bool Decompose();
    void MakeL (BasicMatrix &L);
    void MakeU (BasicMatrix &U);
    bool Solve();
};

class Cholesky : public LinearEquation
{
  public:
    Cholesky();
    BasicMatrix L, LT;
    bool Decompose();
    bool Solve();
};
#endif

