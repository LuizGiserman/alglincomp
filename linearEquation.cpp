/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "linearEquation.h"

LinearEquation::LinearEquation()
{
  unsigned int numVar;
  string str;

  cout << "Enter the number of variables in your system: ";
  getline (cin, str);
  if (!VerifyInput(str, false))
  {
    cout << "Number of variables has to be of type unsigned int\n";
    exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> numVar;	

  this->matrixA.m = numVar;
  cout << "Regarding Matrix A:\n";
  this->matrixA.SetMatrix();
  
  this->matrixB.m = numVar;
  this->matrixB.n = 1;
  cout << "Regarding Matrix B:\n";
  this->matrixB.SetMatrix();

  this->numberVariables = numVar;
  solution.resize(this->numberVariables, 0);

}

bool LinearEquation::ForwardSubstitution(BasicMatrix triangular, BasicMatrix b, BasicMatrix &result)
{
  unsigned int i, j;
  double auxiliar = 0;
  result.m = b.m;
  result.n = 1;
  result.Allocate();
  result.matrix[0][0] = b.matrix[0][0]/triangular.matrix[0][0];

  for (i = 1; i < triangular.m; i++)
  {
    for (j=0; j < i; j++)
    {
      auxiliar += triangular.matrix[i][j] * result.matrix[j][0];
    }
    if (!triangular.matrix[i][i])
    {
      return false;
    }

    result.matrix[i][0] = (b.matrix[i][0] - auxiliar) / triangular.matrix[i][i];
    auxiliar = 0;
  }

  return true;
}

bool LinearEquation::BackwardsSubstitution(BasicMatrix triangular, BasicMatrix b, BasicMatrix &result)
{
  unsigned int i, j;
  double auxiliar = 0;
  result.m = b.m;
  result.n = 1;
  result.Allocate();
  result.matrix[0][0] = b.matrix[0][0]/triangular.matrix[0][0];

  for (i = triangular.m-1; i >= 0; i++)
  {
    for (j=i+1; j < triangular.m-1; j++)
    {
      auxiliar += triangular.matrix[i][j] * result.matrix[j][0];
    }
    if (!triangular.matrix[i][j])
    {
      return false;
    }

    result.matrix[i][0] = (b.matrix[i][0] - auxiliar) / triangular.matrix[i][i];
    auxiliar = 0;
  }

  return true;
}


bool LinearEquation::PrintSolution ()
{
  unsigned int i;

  if(this->solution.size() != this->numberVariables)
  {
    return false;
  }
  
  for (i = 0; i < this->numberVariables; i++)
  {
    cout << "x[" << i+1 << "] = " << this->solution[i] << ";  ";
  }
  cout << endl;
  return true;
}



LU::LU()
{
  this->matrixLU.Copy(this->matrixA);
  this->Decompose();
}

bool LU::Decompose() 
{
  unsigned int k, i, j;

  for (k=0; k < this->matrixLU.n -1; k++)
  {
    for (i=k+1; i < matrixLU.n; i++)
    {
      if (!this->matrixLU.matrix[k][k])
      {
        return false;
      }
      this->matrixLU.matrix[i][k] = this->matrixLU.matrix[i][k]/ this->matrixLU.matrix[k][k];
    }
    for (j = k+1; j < this->matrixLU.n; j++)
    {
      for (i = k+1; i < this->matrixLU.n; i++)
      {
        this->matrixLU.matrix[i][j] = this->matrixLU.matrix[i][j] - this->matrixLU.matrix[i][k] * this->matrixLU.matrix[k][j];
      }
    }
  }
  return true;
}

bool LU::Solve()
{
  BasicMatrix L, U;
  BasicMatrix forwardResult, backwardsResult;

  MakeL(L);
  MakeU(U); 

  L.PrintMatrix();
  U.PrintMatrix();
  if(!this->ForwardSubstitution(L, this->matrixB, forwardResult))
  {
    return false;
  }
  forwardResult.PrintMatrix(); 
  if(!this->BackwardsSubstitution(U, forwardResult, backwardsResult))
  {
    return false;
  }
  backwardsResult.PrintMatrix();
  for (unsigned i = 0; i < backwardsResult.m; i++)
  {
    this->solution[i] = backwardsResult.matrix[i][0];
  }
  return true;
}

void LU::MakeL (BasicMatrix &L)
{
  unsigned int i, j;
  L.m = this-> matrixLU.m;
  L.n = this-> matrixLU.n;
  L.Allocate();

  for (i=0; i < L.m; i++)
  {
    for (j=0; j < L.n; j++)
    {
      if (i==j)
      {
        L.matrix[i][j] = 1;
      }
      else if( i > j)
      {
        L.matrix[i][j] = this->matrixLU.matrix[i][j];
      }
    }
  }
}

void LU::MakeU (BasicMatrix &U)
{
  
  unsigned int i, j;
  U.m = this-> matrixLU.m;
  U.n = this-> matrixLU.n;
  U.Allocate();

  for (i=0; i < U.m; i++)
  {
    for (j=0; j < U.n; j++)
    {
      if (i <= j)
      {
        U.matrix[i][j] = this->matrixLU.matrix[i][j];
      }
    }
  }
}
