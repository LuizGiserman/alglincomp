/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>

using namespace std;

class Utilities 
{
  public:
    Utilities() {};

	  bool VerifyInput(string str, bool canBeNegative, bool isDouble=false);

};


#endif

