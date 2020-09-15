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
  BasicMatrix a;
  a.SetFromFile("matrix.txt");
  cout << "Determinant = " << a.Determinant(a) << endl;
  return 0;

}
