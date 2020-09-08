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

