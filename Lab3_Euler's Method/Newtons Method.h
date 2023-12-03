#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Gauss Method.h"

using namespace std;

namespace precisions {
	const double E1 = 1e-9;
	const double E2 = 1e-9;
}

vector<double> udt(vector<double> u, vector<vector<double>> A, vector<double> b);
vector<double> f(vector<double> y_next, vector<double> y, double h, vector<vector<double>> A, vector<double> b);
vector<vector<double>> calculateMatrixJacobi(vector<double> x, vector<double> y, double h, const double M, const double E1, vector<vector<double>> A, vector<double> b);

vector<double> solveNewtonsMethod(vector<double> x, vector<double> y, double h, vector<vector<double>> A, vector<double> b, const double E1, const double E2, const double M);
