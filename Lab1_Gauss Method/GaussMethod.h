#pragma once
#include <iostream>
#include <vector>

using namespace std;

namespace errors {
	const int MATRIX_WITH_PROPORTIONAL_ROWS = 201;
}

template<typename T>
vector<T> operator / (vector<T>& v1, vector<T>& v2);
template<typename T>
vector<T> operator * (vector<vector<T>>& m1, vector<T>& v2);

template<typename T>
ostream& operator << (ostream& os, vector<vector<T>>& A);
template<typename T>
ostream& operator << (ostream& os, vector<T>& A);
template<typename T>
istream& operator >> (istream& is, vector<vector<T>>& A);
template<typename T>
istream& operator >> (istream& is, vector<T>& A);

template<typename T>
bool isVectorWithSameCoordinates(vector<T> v);
template<typename T>
bool isMatrixWithProportionalityRows(vector<vector<T>> A, vector<T> b);

template<typename T>
vector<T> solveGaussMethod(vector<vector<T>> A, vector<T> b);
template<typename T>
void forwardGaussMethod(vector<vector<T>>& A, vector<T>& b);
template<typename T>
void backwardGaussMethod(vector<vector<T>> A, vector<T> b, vector<T>& x);
template<typename T>
vector<T> calculateResidualVector(vector<vector<T>> A, vector<T> b, vector<T> x);
template<typename T>
T calculateNormaOfRV(vector<T> F);
template<typename T>
T calculateCalcError(vector<T> x1, vector<T> x2);

template<typename T>
vector<T> solveFactorizationLDLT(vector<vector<T>> A, vector<T> b);
