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
  ExOne a;
  cout << "f(277.221) = " << a.Function(277.221) << endl;;
  /*a.Solve();
  cout << "Solution: " << a.solution << endl;*/
  return 0;

}
