/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "basicMatrix.h"

BasicMatrix::BasicMatrix (unsigned int n, unsigned int m) 
{
	
	this->m = m;
	this->n = n;

}

bool BasicMatrix::SetMatrix()
{
  vector<vector<double>> matrix;
	string str;
	if (this->m == 0)
	{
    cout << "Enter the first dimension of the matrix\n";
    getline (cin, str);
    if (!VerifyInput(str, false))
    {
      cout << "Dimension must be of type unsigned int\n";
      exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> m;	
  }
  if(this->n == 0)
  {
    cout << "Enter the second dimension of the matrix\n";
    getline (cin, str);
    if (!VerifyInput(str, false))
    {
      cout << "Dimension must be of type unsigned int\n";
      exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> n;	
	}
	
	matrix.resize(this->m);
  for (unsigned int row = 0; row < m; row++)
	{
		for (unsigned int column = 0; column < n; column++)
		{
			matrix[row].resize(this->n);
			cout << "enter element matrix["<< row + 1 << "][" << column + 1 << "]: "; 
			getline (cin, str);
			if (!VerifyInput(str, true, true))
			{
				cout << "Element must be of type double\n";
			  return false;
      }
			stringstream(str) >> matrix[row][column];
		}
	}

  this->matrix = matrix;
  this->m = m;
  this->n = n;

  return true;
}


void BasicMatrix::PrintMatrix()
{
	for (unsigned int row = 0; row < this->m; row++)
	{
		for (unsigned int column = 0; column < this->n; column++)
			{
					cout << this->matrix[row][column] << "\t";	
			}
		cout << endl;
	}
	cout << endl;
}

bool BasicMatrix::Allocate ()
{
  if (this-> m == 0 || this-> n == 0)
  {
    return false;
  }
  this->matrix.resize(this->m);
  for (unsigned int counter = 0; counter < m; counter++)
  {
    this->matrix[counter].resize(this->n, 0);
  }
  return true;
}

bool BasicMatrix::Add(BasicMatrix matrixB)
{
  if (this->m != matrixB.m || this->n != matrixB.n)
  {
    return false;
  }

  for (unsigned int i = 0; i < matrixB.m; i++)
  {
    for (unsigned int j = 0; j < matrixB.n; j++)
    {
      this->matrix[i][j] += matrixB.matrix[i][j];
    }
  }

  return true;
}

bool BasicMatrix::Add(BasicMatrix matrixB, BasicMatrix &result)
{
  if (this->m != matrixB.m || this->n != matrixB.n)
  {
    return false;
  }

  result.Allocate();
  for (unsigned int i = 0; i < matrixB.m; i++)
  {
    for (unsigned int j = 0; j < matrixB.n; j++)
    {
      result.matrix[i][j] = this->matrix[i][j] + matrixB.matrix[i][j];
    }
  }

  return true;

}

bool BasicMatrix::Cross(BasicMatrix &result, BasicMatrix matrixB)
{
  if (this->n != matrixB.m)
  {
    return false;
  }
  result.n = matrixB.n;
  result.m = this->m;
  result.Allocate();
  for (unsigned int i = 0; i < this->m; i++)
  {
    for (unsigned int j = 0; j < matrixB.n; j++)
    {
      for (unsigned int k = 0; k < matrixB.m; k++)
      {
        result.matrix[i][j] += this->matrix[i][k] * matrixB.matrix[k][j];
      }
    }
  }
  return true;
}

bool BasicMatrix::Copy(BasicMatrix a)
{
  unsigned int i, j;
  
  if(a.m == 0 || a.n == 0)
  {
    return false;
  }

  if (this->m == 0 || this->n == 0)
  {
    this->m = a.m;
    this->n = a.n;
    this->Allocate();
  }

  for (i=0; i < a.m; i++)
  {
    for (j=0; j < a.n; j++)
    {
      this->matrix[i][j] = a.matrix[i][j];
    }
  }
  
  return true;
}
