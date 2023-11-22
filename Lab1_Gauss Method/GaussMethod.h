#pragma once
#include <iostream>
#include <vector>

using namespace std;

namespace errors {
	const int MATRIX_WITH_PROPORTIONAL_ROWS = 201;
}

vector<double> operator / (const vector<double>& v1, const vector<double>& v2);
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2);

ostream& operator << (ostream& os, const vector<vector<double>>& A);
ostream& operator << (ostream& os, const vector<double>& A);
istream& operator >> (istream& is, vector<vector<double>>& A);
istream& operator >> (istream& is, vector<double>& A);

bool isVectorWithSameCoordinates(vector<double> v);
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b);

vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b);
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b);
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x);
vector<double> calculateResidualVector(vector<vector<double>> A, vector<double> b, vector<double> x);
double calculateNormaOfRV(vector<double> F);
double calculateCalcError(vector<double> x1, vector<double> x2);
vector<double> solveFactorizationLDLT(vector<vector<double>> A, vector<double> b);
