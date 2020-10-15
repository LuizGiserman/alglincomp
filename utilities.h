/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <math.h>
#include <vector>
#include "basicMatrix.h"

using namespace std; 

class Utilities 
{
  public:
    Utilities() {};

	  static bool VerifyInput(string str, bool canBeNegative, bool isDouble=false);
    static double GetModule(double number);
    static double FunctionOne(double value);
    static double FunctionTwo(double value);
    static double DerivativeOne(double value);
    static double DerivativeTwo(double value);
};


#endif

