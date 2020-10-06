/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <math.h>

using namespace std;

class Utilities 
{
  public:
    Utilities() {};

	  bool VerifyInput(string str, bool canBeNegative, bool isDouble=false);
    double GetModule(double number);
    double FunctionOne(double value);
    double FunctionTwo(double value);
    double DerivativeOne(double value);
    double DerivativeTwo(double value);
};


#endif

