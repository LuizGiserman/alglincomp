/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef LINEAR_EQUATION_H
#define LINEAR_EQUATION_H

#include "basicMatrix.h"

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
    bool Solve();
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
    bool Solve();
};
#endif

