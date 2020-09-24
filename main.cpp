#include "basicMatrix.h"
#include "linearEquation.h"
#include "mmse.h"

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
  BasicMatrix c;
  c.SetFromFile("matrix.txt");
  c.Jacobi(0.001);
  return 0;

}
