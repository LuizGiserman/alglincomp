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
  unsigned int i;
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
  
  if (result.n == 0 || result.m == 0)
  {
    result.n = matrixB.n;
    result.m = this->m;
  }
  if(result.matrix.empty())
  {
    result.Allocate();
  }

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
  int i, j;
  if (result.m == 0 || result.n == 0)
  {
    result.m = this->n;
    result.n = this->m;
  }
  if (result.matrix.empty())
  {
    result.Allocate();
  }
  for (i=0; i < (int) result.m; i++)
  {
    for (j=0; j < (int) result.n; j++)
    {
      result.matrix[i][j] = this->matrix[j][i];
    }
  }
  return true;
}

double BasicMatrix::Determinant(BasicMatrix matrix)
{
  unsigned int k;
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
        if(j == result.n)
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

bool BasicMatrix::IsDiagonallyDominant()
{
  int size = (int) this->m;
  int i, j;
  double sumLines = 0, sumColumns;

  for (i=0; i < size; i++)
  {
    for (j=0; j < size; j++)
    {
      if (j != i)
      {
        sumLines += GetModule(this->matrix[i][j]);
        sumColumns += GetModule(this->matrix[j][i]);
      }
    }
  
    if (GetModule(this->matrix[i][i]) < sumLines || GetModule(this->matrix[i][i]) < sumColumns)
    {
      return false;
    }
    
    sumLines = 0;
    sumColumns = 0;
  }

  return true;

}

double BasicMatrix::Residue(BasicMatrix auxiliar, BasicMatrix x0)
{
  double normAuxiliar =  auxiliar.VectorNorm();
  double subtraction;

  if (normAuxiliar == -1)
  {
    return -1;
  }

  auxiliar.Subtract(x0);
  subtraction = auxiliar.VectorNorm();
  return subtraction/normAuxiliar;
}

bool BasicMatrix::SetFromFile(string fileName)
{
  
  ifstream myFile;
  vector<double> aux;
  int i, j;
  
  myFile.open (fileName);
  if (!myFile && !(myFile.is_open()))
  {
    std::cout << "Unable to open file" << std::endl;
    exit (ERROR_READING_FILE);
  }

  myFile >> this->m;
  myFile >> this->n;

  if (this->m <= 0 || this->n <= 0)
  {
    return false;
  }

  this->Allocate();
  
  for (i=0; i < (int) this->m; i++)
  {
    for(j=0; j < (int) this->n; j++)
    {
      myFile >> this->matrix[i][j];
    }
  } 
  return true;
}

double BasicMatrix::PowerMethod (double threshold, BasicMatrix &eigenVector)
{
  BasicMatrix x0;
  double oldEigenValue= 1.0;
  double eigenValue;
  double residue;
  x0.m = this->m;
  x0.n = 1;
  int iteration = 0;
  x0.Allocate();
  x0.Fill(1);

  do{
  
    this->Cross(eigenVector, x0);
    eigenValue = eigenVector.matrix[0][0];

    for (auto &rows : eigenVector.matrix)
    {
      for (auto &columns :rows)
      {
        columns /= eigenValue;
      }
    }

    residue = GetModule(eigenValue - oldEigenValue) / GetModule (eigenValue);
    cout << "iteration: " << iteration << endl;
    cout << "old eigenValue = " << oldEigenValue << endl;
    cout << "new eigenValue = " << eigenValue << endl; 
    cout << "old eigenVector:" << endl;
    x0.PrintMatrix();
    cout << "new eigenVector:" << endl;
    eigenVector.PrintMatrix();
    cout << "residue = " << residue << endl;
   
    oldEigenValue = eigenValue;

    x0 = eigenVector;
    eigenVector.Fill(0);    
    iteration++;
  }while(residue > threshold);  

  return eigenValue;  
}

void BasicMatrix::Jacobi (double tolerance)
{
  BasicMatrix P, Pt;
  pair<unsigned int, unsigned int> indicesElement;
  BasicMatrix aNew, aOld, xNew, xOld, auxiliar;
  int iteration = 0;

  aOld.Copy(*this);
  aNew = aOld; 
  xOld.MakeIdentity(this->m);
  xNew = xOld;
  if (!this->IsSymmetric())
  {
    cout << "For Jacobi's iterative method, the matrix has to be symmetric\n";
    exit (ERROR_BAD_INPUT);
  }
  
  while(aNew.VerifyToleranceJacobi(tolerance))
  {
    /*reset values*/
    aOld = aNew;
    xOld = xNew;
    aNew.Clear();
    xNew.Clear();
    auxiliar.Clear();
    
    /*get indices of highest module*/
    indicesElement = aOld.GetIndices();

    /*make matrix P*/
    aOld.MakeP(indicesElement, P);
    
    cout << "iteration: " << iteration << endl;
    cout << "indices i, j = " << indicesElement.first + 1 << ", " << indicesElement.second + 1 << endl;
    cout << "P:\n";
    P.PrintMatrix();
    
    /*Make Ak+1*/
    P.Transpose(Pt);
    /*A(k+1) = Pkt x Ak x Pk*/
    Pt.Cross(auxiliar, aOld);
    auxiliar.Cross(aNew, P);
  
    cout << "aNew: \n";
    aNew.PrintMatrix(); 
    
    /*Make X+1*/
    xOld.Cross(xNew, P);
    cout << "xNew: \n";
    xNew.PrintMatrix();

    iteration++;
  }
  cout << "eigenValues: " << endl;
  xNew.PrintMatrix();
  cout << "eigenVectors: "<< endl;
  aNew.PrintMatrix();
}

bool BasicMatrix::VerifyToleranceJacobi (double tolerance)
{
  int i, j;

  for (i=0; i < (int) this->m; i++)
  {
    for (j=0; j < (int) this->n; j++)
    {
      if (i!=j)
      {
        if (GetModule(this->matrix[i][j]) > tolerance)
        {
          return true;
        }
      }
    }
  }
  return false;
}

pair<unsigned int, unsigned int> BasicMatrix::GetIndices ()
{
  int i, j;
  pair<unsigned int, unsigned int> indicesElement =  make_pair(0, 0);
  double maxValue = 0.0;

  for (i=0; i < (int) this->m; i++)
  {
    for (j=0; j < (int) this->n; j++)
    {
      if (i!=j)
      {
        if (GetModule(this->matrix[i][j]) > maxValue)
        {
          maxValue = GetModule(this->matrix[i][j]);
          indicesElement = make_pair((unsigned int) i, (unsigned int) j);
        }
      }
    }
  }

  return indicesElement;
}

bool BasicMatrix::MakeIdentity(unsigned int size)
{
  int i, j;
  if(size != 0)
  {
    this->m = size;
    this->n = size;
  }

  if (this->m != this->n)
  {
    return false;
  }
  if (this->matrix.empty())
  {
    this->Allocate();
  }

  for(i = 0; i < (int) this->m; i++)
  {
    for (j = 0; j < (int) this->n; j++)
    {
      if (j==i)
      {
        this->matrix[i][j] = 1;
      }
      else
      {
        this->matrix[i][j] = 0;
      }
    }
  }
  return true;
}

bool BasicMatrix::MakeP (pair<unsigned int, unsigned int> indicesElement, BasicMatrix &result)
{
  BasicMatrix identity;
  double phi;
  unsigned int i, j;
  i = indicesElement.first;
  j = indicesElement.second;
  
  if (this->matrix[i][i] != this->matrix[j][j])
  {
    phi = atan((2*this->matrix[i][j])/(this->matrix[i][i] - this->matrix[j][j]))/2;
  }
  else
  {
    phi = PI/4.0;
  }
  cout << "Phi = " << phi << endl; 
  identity.MakeIdentity(this->m);
  identity.matrix[i][i] = cos (phi);
  identity.matrix[i][j] = -1.0*sin(phi);
  identity.matrix[j][j] = cos(phi);
  identity.matrix[j][i] = sin(phi);

  result = identity;
  return true;
}
