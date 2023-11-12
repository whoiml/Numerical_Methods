#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const double E1 = 1e-9;
const double E2 = 1e-9;

namespace errors {
	const int MATRIX_WITH_PROPORTIONAL_ROWS = 201;
}

double f1(double x1, double x2);
double f1dx1(double x1, double x2);
double f1dx2(double x1, double x2);
double f2(double x1, double x2);
double f2dx1(double x1, double x2);
double f2dx2(double x1, double x2);

vector<vector<double>> calculateMatrixJ(double x1, double x2, double k);
void solveNewtonsMethod(double x1, double x2, const double E1, const double E2, const int max_iter, const double M = NULL);