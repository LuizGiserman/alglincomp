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

  cout << "Enter the value A" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "A has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this-> a;

  cout << "Enter the value B" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "B has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this-> b;
  
  cout << "Enter the tolerance" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "Tolerance has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->tol;

}

Bissection::Bissection(double tolerance)
{
  this->tol = tolerance;
}

Bissection::Bissection(double a, double b)
{
  string str; 
  
  cout << "Enter the tolerance" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "Tolerance has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->tol;

}

Bissection::Bissection(double a, double b, double tolerance)
{
  this->tol = tolerance;
  Bissection(a, b);
}


bool Bissection::CheckConditions ()
{
  double valueA, valueB, aux;
  
  valueA = this->Function(this->a);
  valueB = this->Function(this->b);

  cout << "a = " << a << endl << "b = " << b << endl << "valueA = " << valueA << endl << "valueB = " << valueB << endl;
  if (valueA == 0)
  {
    this->solution = this->a;
  }
  else if(valueB == 0)
  {
    this->solution = this->b;
  }
  else if(valueA < 0 && valueB > 0)
  {
  }
  else if(valueB < 0 && valueA > 0)
  {
    aux = this->a;
    this->a = this->b;
    this->b = aux;
  }
  else
  {
    cout <<"One negative and one positive value is needed" << endl;
    exit(ERROR_BAD_INPUT);
  }

  return true;
}

void Bissection::Solve()
{
  
  double a, b, tol, x=0;
  
  a = this-> a;
  b = this-> b;
  tol = this-> tol;
  cout << "a, b, tol : " << a << ", " << b << ", " << tol << endl;
  while (GetModule((b-a)) > tol)
  {
    x = (a+b)/2.0;
    cout << "Module -> " << GetModule(b-a) << ", x = " << x << endl;
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


Newton::Newton()
{  
  string str; 
  
  cout << "Enter the first X" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "X has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->x0;
  

  cout << "Enter the number of iterations" << endl;
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
   cout << "Number of iterations has to be of type int\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->niter;

  cout << "Enter the tolerance" << endl;
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
   cout << "Tolerance has to be of type double\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->tol;

}

Newton::Newton (double x0, double tol)
{
  string str;

  cout << "Enter the number of iterations" << endl;
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
   cout << "Number of iterations has to be of type int\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->niter;

  this->x0 = x0;
  this->tol = tol;

}

Newton::Newton (double x0, int niter, double tol)
{
  this->x0 = x0;
  this->niter = niter;
  this->tol = tol;
}

void Newton::Solve()
{
  int k;
  double newX, oldX;
  oldX = this->x0;
  
  for (k = 0; k < this->niter; k++)
  {
    newX =  oldX - (this->Function(oldX)/this->Derivative(oldX));
    cout << "newX = " << newX << " oldX = " << oldX << endl;
    if (GetModule(newX - oldX) < this->tol)
    {
      cout << "result: "  << newX << endl;
      return;
    }
    oldX = newX;
  }
  cout << "convergence not reached" << endl;
  return;
}

void Newton::SecantSolve()
{
  double deltaX = 0.001;
  double fa, fi;
  double newX, currentX;
  double oldX = this->x0;
  int k;

  currentX = oldX + deltaX;
  fa = this->Function(oldX);
  
  for (k = 0; k < this->niter; k++)
  {
    fi = this->Function(currentX);
    newX =  currentX - fi * (currentX-oldX)/ (fi-fa);
    if(GetModule(newX - currentX) < this->tol)
    {
      cout << "Result : " << currentX << endl;
      return;
    }
    fa = fi;
    oldX = currentX;
    currentX = newX;
  }
  cout << "convergence not reached" << endl;
}


