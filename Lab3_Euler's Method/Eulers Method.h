#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Newtons Method.h"

using namespace std;

vector<double> solveExplicitEulerMethod(vector<double> u, pair<double, double> t_boundaries, double omega, bool showSteps = false);
vector<double> solveImplicitEulerMethod(vector<double> u, pair<double, double> t_boundaries, vector<double> lambda, bool showSteps = false);

// Functions for Explicit Euler method
double u0dt(vector<double> u, double t);
double u1dt(vector<double> u, double t, double x);
vector<double> calculateVectorF(vector<double> y, double t, double x);
vector<double> calculateVectorH(vector<double> f, double h_max, double E);
vector<double> calculateVectorY(vector<double> y, vector<double> f, double h);

// Functions for Implicit Euler method
vector<vector<double>> fillMatrixA(vector<double> lambda);
vector<double> fillVectorB(vector<double> lambda);
vector<double> calculateLocalError(vector<double> y_prev, vector<double> y, vector<double> y_next, double h, double h_prev);
bool isLocalErrorExceedsPermissibleOne(vector<double> E, const vector<double> E_acceptable);
double calculateNextStepOfIntegration(vector<double> E, const vector<double> E_acceptable, double h);