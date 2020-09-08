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
  LU teste;
  if (!teste.Solve())
  {
    cout << "deu false\n";
    return 1;
  }
  cout << "answer:\n";
  teste.PrintSolution();
  return 0;

}
