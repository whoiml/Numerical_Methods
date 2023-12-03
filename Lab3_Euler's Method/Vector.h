#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

vector<double> operator + (const vector<double>& v1, const vector<double>& v2);
vector<double> operator - (const vector<double>& v);
vector<double> operator - (const vector<double>& v1, const vector<double>& v2);
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2);
vector<double> operator * (const vector<double>& v, double scalar);
vector<double> operator / (const vector<double>& v1, const vector<double>& v2);
vector<vector<double>> operator * (const vector<vector<double>>& m, double scalar);

ostream& operator << (ostream& os, const vector<vector<double>>& A);
ostream& operator << (ostream& os, const vector<double>& A);

vector<double> abs(vector<double> v);
vector<double> sqrt(vector<double> v);
double max(vector<double> v);
double min(vector<double> v);