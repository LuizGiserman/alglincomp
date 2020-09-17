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
#include <math.h>
#include <fstream>
#include "utilities.h"

#define ERROR_BAD_INPUT				1
#define ERROR_READING_FILE    2
#define PI 3.141592653589793238462643383279502884197169399

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
  bool Subtract(BasicMatrix matrixB);
  bool Subtract(BasicMatrix matrixB, BasicMatrix &result);
  bool Cross(BasicMatrix &result, BasicMatrix matrixB);
	double Determinant(BasicMatrix m);
  bool Copy(BasicMatrix a);
  bool Transpose(BasicMatrix &result);
  bool Clear();
  void Reduce(unsigned int index, BasicMatrix &result);
  bool IsSymmetric();
  double VectorNorm();
  bool Fill(double value);
  bool IsDiagonallyDominant();
  double Residue(BasicMatrix auxiliar, BasicMatrix x0);
  bool SetFromFile (string fileName);
  double PowerMethod (double treshold, BasicMatrix &eigenVector);
  bool MakeIdentity(unsigned int size);
  void Jacobi(double tolerance);
  private:
    pair<unsigned int, unsigned int> GetIndices ();
    bool VerifyToleranceJacobi(double tolerance);
    bool MakeP (pair <unsigned int, unsigned int> indicesElement, BasicMatrix &result);
};

#endif
