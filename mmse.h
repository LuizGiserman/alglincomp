/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef _MMSE_ 
#define _MMSE_

#include <math.h>
#include <string>
#include "basicMatrix.h"
#include "linearEquation.h"

#define ERROR_NO_SOLUTION     2

class MMSE : public Utilities
{
  public:
    MMSE();
    void Solve();
    unsigned int numVar = 0;
    BasicMatrix X;
    BasicMatrix Y;
    BasicMatrix P, Pt; 
    LU eq;
};

#endif  

