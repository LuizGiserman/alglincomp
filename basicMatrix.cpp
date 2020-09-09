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

bool BasicMatrix::Subtract(BasicMatrix matrixB)
{
  if (this->m != matrixB.m || this->n != matrixB.n)
  {
    return false;
  }

  for (unsigned int i = 0; i < matrixB.m; i++)
  {
    for (unsigned int j = 0; j < matrixB.n; j++)
    {
      this->matrix[i][j] -= matrixB.matrix[i][j];
    }
  }

  return true;
}


bool BasicMatrix::Clear()
{
  unsigned int i, j;
  if (this->matrix.empty())
  {
    return false;
  }

  for (i = 0; i < this->m; i++)
  {
      this->matrix[i].clear();
  }
  this->matrix.clear();
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

bool BasicMatrix::Subtract(BasicMatrix matrixB, BasicMatrix &result)
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
      result.matrix[i][j] = this->matrix[i][j] - matrixB.matrix[i][j];
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

bool BasicMatrix::Transpose(BasicMatrix &result)
{
  unsigned int i, j;
  result.m = this->m;
  result.n = this->n;
  result.Allocate();

  for (i=0; i < result.m; i++)
  {
    for (j=0; j < result.n; j++)
    {
      result.matrix[i][j] = this->matrix[j][i];
    }
  }

  return true;
}

double BasicMatrix::Determinant(BasicMatrix matrix)
{
  unsigned int i, j, k;
  BasicMatrix newMatrix;
  double result = 0.0;
  
  if (matrix.m == 1)
  {
    return matrix.matrix[0][0];
  }

  for (k=0; k < matrix.m; k++)
  {
    matrix.Reduce(k, newMatrix);
    result += matrix.matrix[k][0] * pow(-1, k) * this->Determinant(newMatrix);
  }

  return result;

}

void BasicMatrix::Reduce(unsigned int index, BasicMatrix &result)
{
  
  unsigned int i = 0, j =0;
  unsigned int row, col;
  
  result.Clear();
  result.m = this->m-1;
  result.n = this->n-1;
  result.Allocate();
  
  for (row=0; row < this->m; row++)
  {
    for (col=0; col < this->n; col++)
    {

      if (row != index && col != 0)
      {
        result.matrix[i][j] = this->matrix[row][col];
        j++;
        if(j = result.n-1)
        {
          j = 0;
          i++;
        }
      }

    }
  }
  

}

bool BasicMatrix::IsSymmetric()
{
  unsigned int i, j;
  if (this->m != this-> n)
  {
    return false;
  }

  for (i=0; i < this->m; i++)
  {
    for (j=0; j < this->n; j++)
    {
      if(this->matrix[i][j] != this->matrix[j][i])
      {
        return false;
      }
    }
  }

  return true;
}

double BasicMatrix::VectorNorm()
{
  double auxiliar = 0;
  unsigned int i;
  if(! (this->m == 1 || this->n ==1))
  {
    return -1;
  }

  /*considering normal vector where n=1*/

  for (i=0; i < this->m; i++)
  {
    auxiliar += pow(this->matrix[i][0], 2);
  }
  auxiliar = sqrt(auxiliar);
  
  return auxiliar;
}

bool BasicMatrix::Fill(double value)
{
  unsigned int i, j;

  if(this->m == 0 || this->n == 0)
  {
    return false;
  }

  for (i=0; i < this->m; i++)
  {
    for (j=0; j < this->n; j++)
    {
      this->matrix[i][j] = value;
    }
  }
  return true;
}
