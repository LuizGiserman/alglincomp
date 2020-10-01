/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "utilities.h"
#include "nonLinearSolutions.h"

Bissection::Bissection() 
{
  string str;
  double a, b, tolerance;

  cout << "Enter the value A" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "A has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> a;

  cout << "Enter the value B" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "B has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> b;
  
  if (!CheckConditions(a, b))
  {
    cout <<"One negative and one positive value is needed" << endl;
    exit(ERROR_BAD_INPUT);
  }
  
  cout << "Enter the tolerance" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "Tolerance has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> tolerance;
  this->tol = tolerance;

}

Bissection::Bissection(double tolerance)
{
  this->tol = tolerance;
}

Bissection::Bissection(double a, double b)
{
  if (!CheckConditions(a, b))
  {
    cout <<"One negative and one positive value is needed" << endl;
    exit(ERROR_BAD_INPUT);
  }
}

Bissection::Bissection(double a, double b, double tolerance)
{
  this->tol = tolerance;
  Bissection(a, b);
}


bool Bissection::CheckConditions (double a, double b)
{
  double valueA, valueB;
  
  valueA = Function(a);
  valueB = Function(b);

  if (valueA == 0)
  {
    this->solution = a;
  }
  else if(valueB == 0)
  {
    this->solution = b;
  }
  else if(valueA < 0 && valueB > 0)
  {
    this->a = a;
    this->b = b;
  }
  else if(valueB < 0 && valueA > 0)
  {
    this->a = b;
    this->b = a;
  }
  else
  {
    return false;
  }

  return true;
}

void Bissection::Solve()
{
  
  double a, b, tol, x;
  
  a = this-> a;
  b = this-> b;
  tol = this-> tol;

  while (GetModule(b-a) > tol)
  {
    x = (a+b)/2.0;
    if (Function(x) > 0.0)
    {
      b = x;
    }
    else
    {
      a = x;
    }
  } 

  this->solution = x;
}

double ExOne::Function (double value)
{
  double g = 9.806;
  double k = 0.00341;
  return (double) log(cosh(sqrt(g*k))) - 50.0;
}
