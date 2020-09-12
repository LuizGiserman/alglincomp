#include "basicMatrix.h"
#include "linearEquation.h"

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
  GaussSeidel j;
  j.Solve();
  return 0;

}
