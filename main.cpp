#include "basicMatrix.h"
#include "linearEquation.h"
#include "mmse.h"
#include "nonLinearSolutions.h"


double function1 (vector<double> x)
{
  return x[0] + 2*x[1] - 2;
}

double function2 (vector<double> x)
{
  return x[0]*x[0] + 4 * x[1] * x[1] -4;  
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

  

}
