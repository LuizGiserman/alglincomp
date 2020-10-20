/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "linearEquation.h"

LinearEquation::LinearEquation(string s)
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
  
  if (s.compare("") == 0)
  {
    this->matrixA.m = numVar;
    cout << "Regarding Matrix A:\n";
    cout << "string = " << s << endl;
    this->matrixA.SetMatrix();
     
  }
  else
  {
    this->matrixA.SetFromFile(s);
    if(this->matrixA.m != numVar)
    {
      cout << "Number of variables doesn't match the dimensions of matrix A\n";
      exit (ERROR_BAD_INPUT);
    }
  } 
  
  
  if (this->matrixA.Determinant(this->matrixA) == 0)
  {
    cout << "No consistent solution to this problem\n";
    exit(ERROR_NO_SOLUTION);
  }
  
  this->matrixB.m = numVar;
  this->matrixB.n = 1;
  cout << "Regarding Matrix B:\n";
  this->matrixB.SetMatrix();

  this->numberVariables = numVar;
  solution.resize(this->numberVariables, 0);

}

LinearEquation::LinearEquation(bool empty)
{
  solution.resize(this->numberVariables, 0);
}

LinearEquation::LinearEquation(BasicMatrix a, BasicMatrix b)
{
  this->matrixA = a;
  this->matrixB = b;
  this->numberVariables = b.m;

  if (this->matrixA.Determinant(this->matrixA) == 0)
  {
    cout << "No consistent solution to this problem\n";
    exit(ERROR_NO_SOLUTION);
  }

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

  for (i = (int) triangular.m-1; i >= 0; i--)
  {
    for (j=i+1; j < (int) triangular.m; j++)
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



LU::LU(string s): LinearEquation(s)
{
  if (this->matrixA.m != this->matrixA.n)
  {
    cout << "Matrix must be squared for LU Decomposition" << endl;
    exit (ERROR_BAD_INPUT);
  }
  if (this->matrixA.Determinant(this->matrixA) == 0)
  {
    cout << "Matrix can't be singular for LU Decomposition" << endl;
    exit (ERROR_BAD_INPUT);
  }
  this->matrixLU.Copy(this->matrixA);
  this->Decompose();
}

void LU::Check()
{
  if (this->matrixA.m != this->matrixA.n)
  {
    cout << "Matrix must be squared for LU Decomposition" << endl;
    exit (ERROR_BAD_INPUT);
  }
  if (this->matrixA.Determinant(this->matrixA) == 0)
  {
    cout << "Matrix can't be singular for LU Decomposition" << endl;
    exit (ERROR_BAD_INPUT);
  }
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

  cout << "L: \n";
  L.PrintMatrix();
  cout << "U: \n";
  U.PrintMatrix();
  if(!this->ForwardSubstitution(L, this->matrixB, forwardResult))
  {
    return false;
  }
  if(!this->BackwardsSubstitution(U, forwardResult, backwardsResult))
  {
    return false;
  }
  for (unsigned i = 0; i < backwardsResult.m; i++)
  {
    this->solution[i] = backwardsResult.matrix[i][0];
    cout << "x[" << i << "] = " << solution[i] << " ";
  }

  this->solutionMatrix = backwardsResult;
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

Cholesky::Cholesky(string s): LinearEquation(s)
{
  if (!matrixA.IsSymmetric())
  {
    cout << "Matrix must be symetric for cholesky's decomposition\n";
    exit(ERROR_BAD_INPUT);
  }
  if (!matrixA.IsPositiveDefinite())
  {
    cout << "Matrix must be positive definite\n";
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
  
  BasicMatrix forwardResult, backwardsResult;
  
  if(!this->ForwardSubstitution(L, this->matrixB, forwardResult))
  {
    return false;
  }
  if(!this->BackwardsSubstitution(LT, forwardResult, backwardsResult))
  {
    return false;
  }
  for (unsigned i = 0; i < backwardsResult.m; i++)
  {
    this->solution[i] = backwardsResult.matrix[i][0];
    cout << "x[" << i << "] = " << solution[i] << " ";
  }
  
  return true;
}


Jacobi::Jacobi(string s) : LinearEquation(s)
{
  string str;
  if (!this->matrixA.IsDiagonallyDominant())
  {
    cout << "Matrix A has to be diagonally dominant for Jacobi's method\n" << endl;
    exit(ERROR_BAD_INPUT);
  }
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
  int i, j;
  int size = (int) this->matrixA.m;
  BasicMatrix x0, auxiliar;
  double sum = 0;
  double residue;
  int iteration = 0;

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
        if (j != i)
        {
          sum += this->matrixA.matrix[i][j] * x0.matrix[j][0];
        }
      }
      if (this->matrixA.matrix[i][i] == 0)
      {
        return false;
      }
      auxiliar.matrix[i][0] = (this->matrixB.matrix[i][0] - sum)/this->matrixA.matrix[i][i];
      sum = 0;
    }
    
    residue = this->matrixA.Residue(auxiliar, x0);
    cout << "iteration: " << iteration << endl;
    cout << "residue = " << residue << endl;
    x0 = auxiliar;
    x0.PrintMatrix();
    iteration++;
  }
  while (residue > this->threshold);

  for (i = 0; i < size; i++)
  {
    this->solution[i] = x0.matrix[i][0];
  }

 /* x0.PrintMatrix();*/

  return true;
}

GaussSeidel::GaussSeidel(string s) : LinearEquation(s)
{
  string str;
  if (!this->matrixA.IsDiagonallyDominant())
  {
    if (!this->matrixA.IsSymmetric() || !this->matrixA.IsPositiveDefinite())
    {
      cout << "Matrix A has to be diagonally dominant (or deffinite symmetric) for Gauss-Seidel's method\n" << endl;
      exit(ERROR_BAD_INPUT);
    }
  }
  cout << "Enter the Residual threshold: ";
  getline (cin, str);
  if (!VerifyInput(str, true, true))
  {
    cout << "Threshold has to be of type double\n";
    exit(ERROR_BAD_INPUT);
  }
  stringstream(str) >> this->threshold;

}

bool GaussSeidel::Solve()
{
  
  int i, j;
  int iteration = 0;
  int size = (int) this->matrixA.m;
  BasicMatrix x0, auxiliar;
  double sumNew = 0.0, sumOld = 0.0;
  double residue;

  x0.m = size;
  x0.n = 1;
  x0.Allocate();
  x0.SetMatrix();

  auxiliar.m = size;
  auxiliar.n = 1;
  auxiliar.Allocate();
  
  do
  {
    for (i=0; i < size; i++)
    {
      for (j=0; j < i; j++)
      {
        sumNew += this->matrixA.matrix[i][j] * auxiliar.matrix[j][0]; 
      }
      
      for (j=i+1; j < size; j++)
      {
        sumOld += this->matrixA.matrix[i][j] * x0.matrix[j][0];
      }

      if (this->matrixA.matrix[i][i] == 0)
      {
        return false;
      }
      auxiliar.matrix[i][0] = (this->matrixB.matrix[i][0] - sumNew - sumOld)/this->matrixA.matrix[i][i];
      sumOld = 0.0;
      sumNew = 0.0;
    }

    residue = this->matrixA.Residue(auxiliar, x0); 

    x0 = auxiliar; 
  /*cout << "iteration: " << iteration << endl;
    cout << "residue = " << residue << endl; 
    x0.PrintMatrix();
  */
    auxiliar.Clear();
    auxiliar.Allocate();

    iteration++;
  }while(residue > this->threshold); 
 
  for (i = 0; i < size; i++)
  {
    this->solution[i] = x0.matrix[i][0];
  }

  return true;
}
