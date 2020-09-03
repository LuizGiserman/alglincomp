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

#define ERROR_BAD_INPUT				1
using namespace std;

class BasicMatrix
{
	public:
	BasicMatrix(unsigned int m=0, unsigned int n=0);
	vector<vector<float>> matrix;
	unsigned int m, n;

	float ComputeDeterminant();
	bool VerifyInput(string str, bool canBeNegative, bool isDouble=false);
	void PrintMatrix();

};

#endif
