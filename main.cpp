#include "basicMatrix.h"
#include "linearEquation.h"
#include "mmse.h"
#include "nonLinearSolutions.h"
int main ()
{

  /*LinearEquation eq(lu);
  eq.matrixA.PrintMatrix();
  eq.matrixB.PrintMatrix(); 
  */
  /*BasicMatrix a, b, result;
  a.SetMatrix();
  b.SetMatrix();
  cout << "A:\n";
  a.PrintMatrix();
  cout << "B:\n";
  b.PrintMatrix();
  a.Cross(result, b);
  result.PrintMatrix(); 
 */
  ExOneNewton a(100, 1000, 0.0001);
  cout << "f(100) : " << a.Function(100) << endl;
  cout << "f'(100): " << a.Derivative(100) << endl;
  a.SecantSolve();
  return 0;

}
