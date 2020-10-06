/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef _NON_LINEAR_SOLUTIONS__ 
#define _NON_LINEAR_SOLUTIONS__ 

#include <math.h>
#include <string>
#include "basicMatrix.h"
#include "linearEquation.h"

#define ERROR_NO_SOLUTION     2

class Bissection : public Utilities
{
  public:
    Bissection(double a, double b, double tolerance);
    Bissection(double a, double b);
    Bissection(double tolerance);
    Bissection();
    
    double a, b, solution, tol;

    bool CheckConditions();
    virtual double Function(double value) {return 0.0;};
    void Solve(); 
    
};

class ExOneBissection : public Bissection
{
  public:
    ExOneBissection() {};
    ExOneBissection(double a, double b, double tolerance) : Bissection(a, b, tolerance) {};
    ExOneBissection(double a, double b) : Bissection (a, b) {};
    ExOneBissection(double tolerance) : Bissection (tolerance) {};
    double Function(double value) {return FunctionOne(value);};
};

class ExTwoBissection : public Bissection
{
  public:
    ExTwoBissection() {};
    ExTwoBissection(double a, double b, double tolerance) : Bissection(a, b, tolerance) {};
    ExTwoBissection(double a, double b) : Bissection (a, b) {};
    ExTwoBissection(double tolerance) : Bissection (tolerance) {};
    double Function(double value) {return FunctionTwo(value);};
};

class Newton : public Utilities
{
  public:
    Newton();
    Newton(double x0, double tol);
    Newton (double x0, int niter, double tol);

    int niter;
    double x0, tol;
  
    virtual double Function (double value) {return 0.0;};
    virtual double Derivative (double value) {return 0.0;};
    void Solve();
    void SecantSolve();
};

class ExOneNewton : public Newton
{
  public:
    ExOneNewton() {};
    ExOneNewton(double x0, double tol) : Newton(x0, tol) {}; 
    ExOneNewton(double x0, int niter, double tol) : Newton (x0, niter, tol) {};
    
    double Function(double value) {return this->FunctionOne(value);};
    double Derivative(double value) {return this->DerivativeOne(value);};
};

class ExTwoNewton : public Newton
{
  public:
    ExTwoNewton() {};
    ExTwoNewton(double x0, double tol) : Newton(x0, tol) {}; 
    ExTwoNewton(double x0, int niter, double tol) : Newton (x0, niter, tol) {};
    
    double Function(double value) {return this->FunctionTwo(value);};
    double Derivative(double value) {return this->DerivativeTwo(value);};
};
#endif  

