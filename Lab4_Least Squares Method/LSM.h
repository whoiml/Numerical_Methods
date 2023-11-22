#pragma once
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double operator * (const vector<double> v1, const vector<double> v2);

double sum(vector <double> v);
vector <double> pow(vector<double> v, int degree);
vector <double> log(vector <double> v);

vector <double> solvePowerFunctionApproximation(vector <double> x, vector <double> y);
vector <double> solveLSMApproximation(vector <double> x, vector <double> y);
double calculateStandarDeviation(vector <double> x, vector <double> y, vector <double> parameters);
