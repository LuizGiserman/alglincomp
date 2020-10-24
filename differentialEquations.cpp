/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#include "differentialEquations.h"

EDO::EDO(double (*function) (vector<double>, vector<double>), double t0, double tf, double x0, double delta)
{
    this->func = function;
    this->t0 = t0;
    this->tf = tf;
    this->x0 = x0;
    this->delta = delta;
}

void EDO::Solve()
{
    int steps, i;
    double tNew, tOld;
    double xNew, xOld;
    
    steps = (int) (tf-t0) / delta;
    
    tOld = this->t0;
    xOld = this->x0;

    for (i = 0; i < steps; i++)
    {
        tNew = (i + 1) * delta;
        xNew = this->Method(xOld, tOld);
        cout << "x : " << xNew << "\t t: " << tNew << endl;
        xOld = xNew;
        tOld = tNew;

    }

    cout << "x : " << xNew << endl << "t : " << tNew << endl;
}

double EDOEuler::Method(double xOld, double tOld)
{
    vector<double> x {xOld};
    vector<double> t {tOld};
    return (xOld + delta * this->func (t, x));
}   

double EDORK2::Method(double xOld, double tOld)
{
    vector<double> x {xOld};
    vector<double> t {tOld};
    double k1, k2;

    k1 = this->func(t, x);
    t[0] += this->delta;
    x[0] = x[0] + this->delta * k1;
    k2 = this->func(t, x);
    return xOld + delta / 2 * (k1 + k2);
}

double EDORK4::Method(double xOld, double tOld)
{
    vector<double> x {xOld};
    vector<double> t {tOld};
    double k1, k2, k3, k4;
    
    k1 = this->func(t, x);
    
    t[0] += this->delta / 2;
    x[0] = x[0] + this->delta / 2 * k1;
    k2 = this->func(t, x);
    
    x[0] = xOld + this->delta / 2 * k2;
    k3 = this->func(t, x);
    
    t[0] += this->delta/ 2;
    x[0] = xOld + delta * k3;
    k4 = this->func(t, x);
    
    return (xOld + this->delta / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
}

EDOSecondOrder::EDOSecondOrder(double (*function) (vector<double>, vector<double>, double), double t0, double tf, double x0, double derivative, double delta)
{
    this->func = function;
    this->t0 = t0;
    this->tf = tf;
    this->x0 = x0;
    this->derivative = derivative;
    this->delta = delta;
}


void EDOSecondOrder::Solve()
{
    int steps, i;
    double aux;

    vector<double> t {this->t0};
    vector<double> x {this->x0};

    steps = (int) (this->tf - this->t0) / this->delta;


    for (i = 0; i < steps; i++)
    {
        t.push_back((i+1) * this->delta);
        this->Method(t, x, &aux, i);
        cout << "x : " << x[i] << "\t" << "t : " << t[i] << endl;
    }
    cout << "x : " << x[i] << "\t" << "t : " << t[i] << endl;

}

void EDO_SO_Taylor::Method(vector<double> &t, vector<double> &x, double *aux, int index)
{
    double temp;
    temp = this->func (vector<double> {t[index]}, vector<double>{x[index]}, *aux);
    x.push_back(x[index] + (*aux * this->delta) + ((temp * this->delta * this->delta) / 2));
    *aux += temp * delta;
}

void EDO_RK_Nystrom::Method(vector<double> &t, vector<double> &x, double *aux, int index)
{
    double k1, k2, k3, k4;
    double q, l;

    k1 = this->delta /2 * func (vector<double> {t[index]}, vector<double> {x[index]}, *aux);
    q = this->delta /2 * (*aux + k1 / 2);
    k2 = this->delta /2 * func(vector<double> {t[index] + this->delta/2}, vector<double> {x[index] + q }, *aux + k1);
    k3 = this->delta /2 * func(vector<double> {t[index] + this->delta/2}, vector<double> {x[index] + q }, *aux + k2);
    l = this->delta * (*aux + k3);
    k4 = this->delta/2 * func(vector<double> {t[index]}, vector<double> {x[index] + l}, *aux + 2*k3);

    x.push_back(x[index] + this->delta * (*aux + (k1 + k2 + k3) /3));
    *aux += (k1 + 2*k2 + 2*k3 + k4) / 3;
}