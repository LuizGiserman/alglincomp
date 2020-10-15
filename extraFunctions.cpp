/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "extraFunctions.h"

BasicMatrix Jacobian(vector <double (*)(vector <double> )> listFunctions, vector<double> firstSolution)
{
    int i, j;
    int numberFunctions = listFunctions.size();
    int numberSolutions = firstSolution.size();
    BasicMatrix result ((unsigned int) numberFunctions, (unsigned int) numberSolutions);
    result.Allocate();
    for (i = 0; i < numberFunctions; i++)
    {
        for (j=0; j < numberSolutions; j++)
        {
            result.matrix[i][j] = PartialDerivative(listFunctions[i], firstSolution, j);
        }
    }

    return result;
}

double PartialDerivative(double (*function)(vector<double>), vector<double> firstSolution, int index)
{
	double numerator, denominator, result;
	double delta = pow(10, -10);
    double aux = (*function)(firstSolution);

    vector <double> newSolution = firstSolution;
    newSolution[index] += delta;

    numerator = (*function)(newSolution) - aux;
    denominator = delta;
    result = numerator/denominator;

    return result;
}

double Derivative1D(double (*function)(double value), double value)
{
	double delta, numerator, denominator;

    delta = pow(10,-10);
    numerator = (*function)(value + delta) - (*function)(value);
    denominator = delta;
    return numerator/denominator;
}