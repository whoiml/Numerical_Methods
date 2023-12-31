#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

namespace errors {
	const int MATRIX_WITH_PROPORTIONAL_ROWS = 201;
}

namespace precisions {
	const double E1 = 1e-9;
	const double E2 = 1e-9;
}

/* f1 and f2 are systems of equations that we solve
	dx1 and dx2 are derivatives of them */
double f1(double x1, double x2);
double f1dx1(double x1, double x2);
double f1dx2(double x1, double x2);
double f2(double x1, double x2);
double f2dx1(double x1, double x2);
double f2dx2(double x1, double x2);

vector<vector<double>> calculateJacobiMatrix(double x1, double x2, double k);
void solveNewtonsMethod(double x1, double x2, const double E1, const double E2, const int max_iter, const double M = NULL);
