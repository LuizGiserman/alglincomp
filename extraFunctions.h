/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef EXTRA_H
#define EXTRA_H

#include <string>
#include <math.h>
#include <vector>
#include "basicMatrix.h"

using namespace std; 


void Jacobian (vector <double (*)(vector <double> )> listFunctions, BasicMatrix firstSolution, BasicMatrix &J);
double PartialDerivative(double (*function)(vector<double>), vector<double> firstSolution, int index);
double Derivative1D(double (*function)(double value), double value);
double DerivativeForward (double (*function)(double value), double x, double delta);
double DerivativeCentral (double (*function)(double value), double x, double delta);
double DerivativeBackwards (double (*function)(double value), double x, double delta);


class Richard
{
    public:
    Richard (double (*function)(double value), double x, double delta, double p);
    double (*function)(double value);
    double x, delta, p;
    virtual double Derivative (){return 0.0;};
    double Compute();
};
#endif

class ForwardRichard : public Richard 
{
    ForwardRichard (double (*function)(double value), double x, double delta, double p) : Richard (function, x, delta, p) {};
    double Derivative(){return DerivativeForward(this->function, this->x, this->delta);};
};

class CentralRichard : public Richard 
{
    CentralRichard (double (*function)(double value), double x, double delta, double p) : Richard (function, x, delta, p) {};
    double Derivative(){return DerivativeCentral(this->function, this->x, this->delta);};
};

class BackwardsRichard : public Richard 
{
    BackwardsRichard (double (*function)(double value), double x, double delta, double p) : Richard (function, x, delta, p) {};
    double Derivative(){return DerivativeBackwards(this->function, this->x, this->delta);};
};