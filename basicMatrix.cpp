/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "basicMatrix.h"

BasicMatrix::BasicMatrix (unsigned int n, unsigned int m) 
{
	
	vector<vector<float>> matrix;
	string str;
	if (m == 0 && n == 0)
	{
    cout << "Enter the first dimension of the matrix\n";
    getline (cin, str);
    if (!VerifyInput(str, false))
    {
      cout << "Dimension must be of type unsigned int\n";
      exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> m;	

    cout << "Enter the second dimension of the matrix\n";
    getline (cin, str);
    if (!VerifyInput(str, false))
    {
      cout << "Dimension must be of type unsigned int\n";
      exit(ERROR_BAD_INPUT);
    }
    stringstream(str) >> n;	
	}

	matrix.resize(m);
	for (unsigned int row = 0; row < m; row++)
	{
		for (unsigned int column = 0; column < n; column++)
		{
			matrix[row].resize(n);
			cout << "enter element a["<< row + 1 << "][" << column + 1 << "]: "; 
			getline (cin, str);
			if (!VerifyInput(str, true, true))
			{
				cout << "Dimension must be of type double\n";
				exit(ERROR_BAD_INPUT);
			}
			stringstream(str) >> matrix[row][column];
		}
	}
	
	this->m = m;
	this->n = n;
	this->matrix = matrix;

}

bool BasicMatrix::VerifyInput(string str, bool canBeNegative, bool isDouble)
{
	for (auto const &x: str)
	{
		if ( x < '0' || x > '9')
		{
			if(!isDouble && !canBeNegative)
			{
				return false;
			}
			if (x != '.')
			{
				if (x != '-')
				{
					return false;
				}
			}
		}
	}

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

