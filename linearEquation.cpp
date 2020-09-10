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
  
  if (this->matrixA.Determinant(this->matrixA) == 0)
  {
    cout << "No consistent solution to this problem\n";
    exit(ERROR_NO_SOLUTION);
  }

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
  int i, j;
  double auxiliar = 0;
  result.m = triangular.m;
  result.n = 1;
  result.Allocate();
  result.matrix[0][0] = b.matrix[0][0]/triangular.matrix[0][0];

  for (i = triangular.m-1; i >= 0; i--)
  {
    for (j=i+1; j < triangular.m; j++)
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

Cholesky::Cholesky()
{
  if (!matrixA.IsSymmetric())
  {
    cout << "Matrix must be symetric for cholesky's decomposition\n";
    exit(ERROR_BAD_INPUT);
  }
  this->L.m = this->matrixA.m;
  this->L.n = this->matrixA.n;
  this->L.Allocate();
  this->Decompose();
}


bool Cholesky::Decompose()
{
  unsigned int i, j, k;
  double auxiliar = 0.0;
  for (i=0; i < this->L.m; i++)
  {
    for (k=0; k < i; k++)
    {
      auxiliar += pow(L.matrix[i][k], 2);
    }
  
    this->L.matrix[i][i] = sqrt(this->matrixA.matrix[i][i] - auxiliar);
    auxiliar = 0;
  
    for (j=i+1; j < this->L.m; j++)
    {
      for (k=0; k < i; k++)
      {
        auxiliar += this->L.matrix[i][k] * this->L.matrix[j][k];
      } 

      L.matrix[j][i] = (this->matrixA.matrix[i][j] - auxiliar) / this->L.matrix[i][i];
      auxiliar = 0.0;
    }
  }
  
  this->L.Transpose(this->LT);
  return true;
}

bool Cholesky::Solve()
{
  return true;

}

Jacobi::Jacobi()
{
  string str;
  cout << "Enter the Residual threshold: ";
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
    cout << "Threshold has to be of type double\n";
    exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->threshold;
}

bool Jacobi::Solve()
{
  unsigned int i, j;
  unsigned int size = this->matrixA.m;
  BasicMatrix x0, auxiliar;
  double sum = 0;
  double residue;

  x0.m = size;
  x0.n = 1;
  x0.Allocate();
  x0.Fill(1);
  auxiliar.m = size;
  auxiliar.n = 1; 
  
  do 
  {
    auxiliar.Clear();
    auxiliar.Allocate();
    
    for (i=0; i < size; i++)
    {
      for (j=0; j < size; j++)
      {
        sum = this->matrixA.matrix[i][j] * x0.matrix[j][0];
      }
      if (this->matrixA.matrix[i][i] == 0)
      {
        return false;
      }
      auxiliar.matrix[i][0] = (this->matrixB.matrix[i][0] - sum)/this->matrixA.matrix[i][i];
      sum = 0;
    }
    residue = this->Residue(auxiliar, x0);
    cout << "residue = " << residue << endl;
    x0 = auxiliar;
    x0.PrintMatrix();
    break;
  }
  while (residue> this->threshold);

  for (i = 0; i < size; i++)
  {
    this->solution[i] = x0.matrix[i][0];
  }

  x0.PrintMatrix();  
  return true;
}

double Jacobi::Residue(BasicMatrix auxiliar, BasicMatrix x0)
{
  double normAuxiliar =  auxiliar.VectorNorm();
  double subtraction;

  if (normAuxiliar == -1)
  {
    return -1;
  }

  auxiliar.Subtract(x0);
  auxiliar.PrintMatrix();
  subtraction = auxiliar.VectorNorm();
  cout << "subtraction =  " << subtraction << endl;
  cout << "normAuxiliar = " << normAuxiliar << endl;
  cout << "division = " << subtraction/normAuxiliar << endl; 
  return subtraction/normAuxiliar;
}
