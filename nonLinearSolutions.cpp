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

InverseInterpolation::InverseInterpolation ()
{
  int i;
  string str; 
  string elements [] = {"first", "second", "third"};
 
  this->xInicial.resize(3, 0);

  for (i=0; i < 3; i++)
  {
    cout << "Enter the " << elements[i] << " X" << endl; 
    getline (cin, str);
    if (!VerifyInput(str, true, true))
    {
      cout << "X has to be of type double\n";
      exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> this->xInicial[i];

  }
  
  cout << "Enter the number of iterations" << endl;
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
   cout << "Number of iteration has to be of type int\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->niter;
  
}

InverseInterpolation::InverseInterpolation (double x1, double x2, double x3)
{
  string str;

  this->xInicial.resize(3, 0);
  this->xInicial[0] = x1;
  this->xInicial[1] = x2;
  this->xInicial[2] = x3;
  sort(this->xInicial.begin(), this->xInicial.end());

  cout << "Enter the number of iterations" << endl;
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
   cout << "Number of iteration has to be of type int\n";
   exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->niter;
  
}

InverseInterpolation::InverseInterpolation (double x1, double x2, double x3, int niter)
{
  this->xInicial.resize(3, 0);
  this->xInicial[0] = x1;
  this->xInicial[1] = x2;
  this->xInicial[2] = x3;
  sort(this->xInicial.begin(), this->xInicial.end());
  this->niter = niter;
}

double InverseInterpolation::LaGrange(vector <double> x)
{
  vector<double> y(3, 0);

  y[0] = this->Function(x[0]);
  y[1] = this->Function(x[1]);
  y[2] = this->Function(x[2]);
 
  return (y[1]*y[2]*x[0]) / ((y[0]-y[1]) * (y[0]-y[2])) +
         (y[0]*y[2]*x[1]) / ((y[1]-y[0]) * (y[1]-y[2])) +
         (y[0]*y[1]*x[2]) / ((y[2]-y[0]) * (y[2]-y[1]));
}

pair<int, double> InverseInterpolation::GetElementAndIndex(vector <double> x) 
{
  pair <int, double> result;
  int count = 0;
  int maxi = 0;
  double max = 0;
  double test;
  for (auto const &aux : x)
  {
    if( (test = GetModule(this->Function(aux))) > max)
    {
      max = test;
      maxi = count;
    }
    count++;
  }

  result = make_pair(maxi, max);
  return result;
}

void InverseInterpolation::Solve()
{
  int k;
  double auxNew, x0 = 10e36;
  vector <double> x;
  pair < int, double > indexValue;
  
  x = this->xInicial;

  for (k = 1; k <= this->niter; k++)
  {
    auxNew = this->LaGrange(x);

    if (GetModule(auxNew - x0) < this->tol)
    {
      cout << "Result : " << auxNew << endl;
      return;
    }
    
    indexValue = GetElementAndIndex(x);
    x[indexValue.first] = auxNew;
    sort(x.begin(), x.end(), [this](double& w1, double& w2){return (this->Function(w1) < this->Function(w2));});
    x0 = auxNew;
    cout << endl;
  } 
   
}


NonLinearEquations::NonLinearEquations(vector<double> firstSolution)
{
  this->firstSolution = firstSolution;
}

NLE_Newton::NLE_Newton(vector <double (*)(vector <double> )> listFunctions, vector<double> firstSolution) : NonLinearEquations(firstSolution) 
{
  this->listFunctions = listFunctions;
}

void NLE_Newton::Solve()
{
  BasicMatrix J;
  int niter = 100;
  int k;
  vector<double> F;
  vector<double> solutionX;
  solutionX = this->firstSolution;
  for (k = 0; k < niter; k++)
  {
    J = Jacobian (this->listFunctions, solutionX);
    F = GetF(solutionX);
    

  }
  
}

vector<double> NLE_Newton::GetF(vector<double> solutionX)
{
  int i = 0;
  vector<double> result (this->listFunctions.size(), 0);
  for (auto const &function: this->listFunctions)
  {
    result[i] = function(solutionX);
  }

  return result;
}