/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef _NON_LINEAR_SOLUTIONS__ 
#define _NON_LINEAR_SOLUTIONS__ 

#include <math.h>
#include <string>
#include <algorithm>
#include "linearEquation.h"
#include "extraFunctions.h"
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
    ExOneBissection() {this->CheckConditions();};
    ExOneBissection(double a, double b, double tolerance) : Bissection(a, b, tolerance) {this->CheckConditions();};
    ExOneBissection(double a, double b) : Bissection (a, b) {this->CheckConditions();};
    ExOneBissection(double tolerance) : Bissection (tolerance) {this->CheckConditions();};
    double Function(double value) {return FunctionOne(value);};
};

class ExTwoBissection : public Bissection
{
  public:
    ExTwoBissection() {this->CheckConditions();};
    ExTwoBissection(double a, double b, double tolerance) : Bissection(a, b, tolerance) {this->CheckConditions();};
    ExTwoBissection(double a, double b) : Bissection (a, b) {this->CheckConditions();};
    ExTwoBissection(double tolerance) : Bissection (tolerance) {this->CheckConditions();};
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

class InverseInterpolation : public Utilities
{
  public:
    InverseInterpolation();
    InverseInterpolation(double x1, double x2, double x3);
    InverseInterpolation(double x1, double x2, double x3, int niter);
    void Solve();

    vector<double> xInicial;
    int niter;
    double tol =0.0001;

    virtual double Function (double value) {return 0.0;};
    double LaGrange(vector <double> x);
    pair<int, double> GetElementAndIndex(vector <double> x);
};

class ExOneII : public InverseInterpolation
{
  public:
    ExOneII() {};
    ExOneII(double x1, double x2, double x3) : InverseInterpolation(x1, x2, x3) {};
    ExOneII(double x1, double x2, double x3, int niter) : InverseInterpolation(x1, x2, x3, niter) {}; 

    double Function(double value) {return this->FunctionOne(value);};
};

class NonLinearEquations : public Utilities
{
  public:
    NonLinearEquations(vector<double> firstSolution);
    vector<double> firstSolution;
};

class NLE_Newton : public NonLinearEquations
{
  public:
  NLE_Newton(vector <double (*)(vector <double> )> listFunctions, vector<double> firstSolution);
  vector <double (*)(vector <double> )> listFunctions;
  vector<double> GetF(vector<double> value);
  void Solve();
};

class NLE_Broyden : public NonLinearEquations
{
  public:

};

