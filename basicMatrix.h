/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef BASIC_MATRIX_H	
#define BASIC_MATRIX_H	

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "utilities.h"

#define ERROR_BAD_INPUT				1
using namespace std;

class BasicMatrix : public Utilities
{
	public:
	BasicMatrix(unsigned int m=0, unsigned int n=0);
	vector<vector<double>> matrix;
	unsigned int m, n;

  bool SetMatrix();
  bool Allocate ();
	void PrintMatrix();
  bool Add(BasicMatrix matrixB);
  bool Add(BasicMatrix matrixB, BasicMatrix &result);
  bool Cross(BasicMatrix &result, BasicMatrix matrixB);
	float ComputeDeterminant();
  bool Copy(BasicMatrix a);
};

#endif
