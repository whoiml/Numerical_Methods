#pragma once
#include "Newtons Method.h"

vector<double> operator / (const vector<double>& v1, const vector<double>& v2);
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2);
vector<double> operator - (const vector<double>& v);
ostream& operator << (ostream& os, const vector<vector<double>>& A);
ostream& operator << (ostream& os, const vector<double>& A);

bool isVectorWithSameCoordinates(vector<double> v);
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b);
vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b);
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b);
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x);