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
  ExOneII a(280, 290, 300, 100);
  a.Solve();
  return 0;

}
