#include "basicMatrix.h"
#include "linearEquation.h"
#include "mmse.h"
#include "nonLinearSolutions.h"
#include "integrals.h"

double function1 (vector<double> x)
{
  return x[0] + 2*x[1] - 2;
}

double function2 (vector<double> x)
{
  return x[0]*x[0] + 4 * x[1] * x[1] -4;  
}

double function3 (vector<double>x)
{
  return exp(-pow(x[0], 2)/2)/sqrt(2 * PI);
}

int main ()
{

  // vector <double> first = {2, 3};

  // vector <double (*)(vector <double> )> funcs;
  // funcs.push_back(function1);
  // funcs.push_back(function2);
  // cout << "pushed back" << endl;
  // // NLE_Newton teste(funcs, first);
  // // cout << "Created Object" << endl;
  // // teste.firstSolution.PrintMatrix();
  // // teste.Solve();
  // NLE_Broyden teste (funcs, first);
  // teste.Solve();
  // return 0;

  IntegralPolinomial teste (function3, -INF, 1, 10);
  cout << "Result : " << teste.Integrate() << endl;
  return 0;

}
