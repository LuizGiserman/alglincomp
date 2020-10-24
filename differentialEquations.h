/*
* Universidade Federal do Rio de Janeiro
* Author: Luiz Giserman
* Class: Numeric Linear Algebra
*/

#ifndef DIFFERENTIAL_EQUATIONS_H
#define DIFFERENTIAL_EQUATIONS_H

#include <vector>
#include <iostream>


using namespace std;

class EDO
{
    public:
    EDO(double (*function) (vector<double>, vector<double>), double t0, double tf, double x0, double delta);
    double (*func) (vector<double>, vector<double>);
    double t0, tf, x0, delta; 
    void Solve();
    virtual double Method (double xOld, double tOld) {return 0.0;};
};

class EDOEuler : public EDO
{
    public:
    EDOEuler(double (*function) (vector<double>, vector<double>), double t0, double tf, double x0, double delta) : EDO (function, t0, tf, x0, delta) {};
    double Method(double xOld, double tOld);

};

class EDORK2 : public EDO
{
    public:
    EDORK2(double (*function) (vector<double>, vector<double>), double t0, double tf, double x0, double delta) : EDO (function, t0, tf, x0, delta) {};
    double Method(double xOld, double tOld);
};

class EDORK4 : public EDO
{
    public:
    EDORK4(double (*function) (vector<double>, vector<double>), double t0, double tf, double x0, double delta) : EDO (function, t0, tf, x0, delta) {};
    double Method(double xOld, double tOld);
};


class EDOSecondOrder
{
    public:
    EDOSecondOrder (double (*function) (vector<double>, vector<double>, double), double t0, double tf, double x0, double derivative, double delta);
    double (*func) (vector<double>, vector<double>, double);
    double t0, tf, x0, derivative, delta;
    void Solve();
    virtual void Method(vector<double> &t, vector<double> &x, double *aux, int index) {return;};
};

class EDO_SO_Taylor : public EDOSecondOrder
{
    public:
    EDO_SO_Taylor (double (*function) (vector<double>, vector<double>, double), double t0, double tf, double x0, double derivative, double delta) : EDOSecondOrder (function, t0, tf, x0, derivative, delta) {};
    void Method(vector<double> &t, vector<double> &x, double *aux, int index) override;
};

class EDO_RK_Nystrom : public EDOSecondOrder
{
    public:
    EDO_RK_Nystrom (double (*function) (vector<double>, vector<double>, double), double t0, double tf, double x0, double derivative, double delta) : EDOSecondOrder (function, t0, tf, x0, derivative, delta) {};
    void Method(vector<double> &t, vector<double> &x, double *aux, int index) override;
};

#endif