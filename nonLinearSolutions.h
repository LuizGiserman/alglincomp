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

    bool CheckConditions(double a, double b);
    virtual double Function(double value) {return 0.0;};
    void Solve(); 
    
};

class ExOne : public Bissection
{
  public:
    double Function(double value);
};

#endif  

