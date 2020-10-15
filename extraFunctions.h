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


BasicMatrix Jacobian (vector <double (*)(vector <double> )> listFunctions, vector<double> firstSolution);
double PartialDerivative(double (*function)(vector<double>), vector<double> firstSolution, int index);
double Derivative1D(double (*function)(double value), double value);


#endif