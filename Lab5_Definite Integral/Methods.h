#pragma once
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

const double PI = 3.1415926535;
const double E4 = 1e-4;
const double E5 = 1e-5;

double f_v19(double x);
double f_v32(int i, int j, int n);

double solveTrapezoidMethod(double a, double b, double(*f)(double), double e);
double solveSimpsonsMethod(double a, double b, double(*f)(double), double E);
double solveSimpsonsCubatureMethod(double a, double A, double b, double B, double(*f)(int, int, int), double e);