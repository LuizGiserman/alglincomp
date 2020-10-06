/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "utilities.h"

bool Utilities::VerifyInput(string str, bool canBeNegative, bool isDouble)
{
	for (auto const &x: str)
	{
		if ( x < '0' || x > '9')
		{
			if(!isDouble && !canBeNegative)
			{
				return false;
			}
			if (x != '.')
			{
				if (x != '-')
				{
					return false;
				}
			}
		}
	}

	return true;
}

double Utilities::GetModule(double number)
{
  if (number < 0)
  {
    return -1 * number;
  }

  return number;
}


double Utilities::FunctionOne (double value)
{
  double g = 9.806;
  double k = 0.00341;
  return (double) log(cosh(value * sqrt(g*k))) - 50.0;
}

double Utilities::DerivativeOne (double value)
{
  double g = 9.806;
  double k = 0.00341;
  double root = sqrt(g*k);
  return (double) root * tanh(root * value);
}

double Utilities::FunctionTwo (double value)
{
  return 4.0*cos(value) - exp(2.0*value);
}

double Utilities::DerivativeTwo (double value)
{
  return -2.0*(exp(2*value) + 2*sin(value));
}
