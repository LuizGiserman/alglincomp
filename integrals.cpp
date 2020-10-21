/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "linearEquation.h"
#include "integrals.h"

Integral::Integral(double (*function) (vector<double>), double a, double b, int numberPoints)
{
    this->function = function;
    this->a = a;
    this->b = b;
    this->numberPoints = numberPoints;
}

double Integral::HermiteFunction (vector<double> value)
{
    return this->function(value) / (exp(-value[0]));
}

double Integral::LaguerreFunction (vector<double> value)
{
    return this->function(value) / (exp(-pow(value[0], 2)));
}

double Integral::Integrate()
{
    vector<double> pointsHermite = HERMITE[this->numberPoints-2].first;
    vector<double> weightHermite = HERMITE[this->numberPoints-2].second;

    vector<double> pointsLaguerre = LAGUERRE[this->numberPoints-2].first;
    vector<double> weightsLaguerre = LAGUERRE[this->numberPoints-2].second;

    double hermiteSum = 0;
    double laguerreSum = 0;

    for (int i=0; i < this->numberPoints; i++)
    {
        hermiteSum += this->HermiteFunction(vector<double>{pointsHermite[i]}) * weightHermite[i];
        laguerreSum += this->LaguerreFunction(vector<double>{pointsLaguerre[i]}) * weightsLaguerre[i];
    }
    cout << "HermiteSum : " << hermiteSum << endl;
    cout << "LaguerreSum : " << laguerreSum << endl;
    
    if (this->a == -INF)
    {
        if (this->b == INF)
        {
            return hermiteSum;
        }
        if (this->b >= 0)
        {   
            this->a = 0;
            return hermiteSum - laguerreSum + Integration();
        }
        this->a = b;
        this->b = 0;
        return hermiteSum - laguerreSum - Integration();
    }
    
    if(this-> b == INF)
    {
        if (this-> a == 0)
        {
            return laguerreSum;
        }
        if (this->a > 0)
        {
            this->b = a;
            this->a = 0;
            return laguerreSum - Integration();
        }

        this->b = 0;
        return laguerreSum + Integration();
    }

    return Integration();

}

double IntegralPolinomial::Integration()
{
    int i, j;
    double deltaX = Utilities::GetModule(this->b - this->a) /(this->numberPoints-1);
    vector<double> vb, incognitos;
    BasicMatrix vandermont, matrixB;
    double result = 0;
    
    vb.reserve(this->numberPoints);
    incognitos.reserve(this->numberPoints);

    for (i = 1; i <= this->numberPoints; i++)
    {
        incognitos.push_back(this->a+(i-1)*deltaX);
        vb.push_back((pow(this->b, i)-pow(this->a, i))/i);
    }

    vandermont.m = this->numberPoints;
    vandermont.n = this->numberPoints;
    vandermont.Allocate();

    for(i = 0; i < (int) vandermont.m; i++)
    {
        for (j = 0; j < (int) vandermont.n; j++)
        {
            vandermont.matrix[i][j] = pow(incognitos[j], i);
        }
    }

    matrixB.m = vb.size();
    matrixB.n = 1;
    matrixB.Allocate();
    i = 0;
    for (auto const &element: vb)
    {
        matrixB.matrix[i][0] = element;
        i++;
    }
    LU vLU(vandermont, matrixB);
    vLU.Solve();
    for (i = 0; i < this->numberPoints; i++)
    {
        result += vLU.solution[i] * this->function(vector<double>{incognitos[i]});
    }
    
    return result;
}

double IntegralQuadr::Integration()
{
    int i;
    vector<double> pointsLegendre, weightsLegendre, incognitos;
    pointsLegendre = LEGENDRE[this->numberPoints-2].first;
    weightsLegendre = LEGENDRE[this->numberPoints-2].second;
    double L = this->b - this->a;
    double integralSum = 0;
    incognitos.reserve(this->numberPoints);
    
    for (i=0; i < this->numberPoints; i++)
    {
        incognitos.push_back( (this->a + this->b + pointsLegendre[i] * L) / 2);
        integralSum += function(vector<double>{incognitos[i]}) * weightsLegendre[i];
    }

    return integralSum * L / 2;
}