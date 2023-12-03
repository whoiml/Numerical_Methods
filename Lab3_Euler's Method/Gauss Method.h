#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Vector.h"

using namespace std;

namespace errors {
	const int MATRIX_WITH_PROPORTIONAL_ROWS = 201;
}

bool isVectorWithSameCoordinates(vector<double> v);
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b);
vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b);
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b);
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x);