#include "LSM.h"
#include "Gauss Method.h"

double operator * (const vector<double> v1, const vector<double> v2) {
	if (v1.size() != v2.size()) {
		throw invalid_argument("Vectors must have the same size.");
	}

	double result = 0.0;
	for (size_t i = 0; i < v1.size(); i++) {
		result += v1[i] * v2[i];
	}

	return result;
}

double sum(vector <double> v) {
	double sum = NULL;
	for (int i = 0; i < v.size(); i++)
		sum += v[i];
	return sum;
}
vector <double> pow(vector<double> v, int degree) {
	vector <double> result(v.size(), 1.0);
	for (int i = 0; i < v.size(); i++) {
		result[i] = pow(v[i], degree);
	}
	return result;
}
vector <double> log(vector <double> v) {
	vector <double> result(v.size());
	for (int i = 0; i < v.size(); i++)
		result[i] = log(v[i]);
	return result;
}

vector <double> solvePowerFunctionApproximation(vector <double> x, vector <double> y) {

	// y=ax^b => ln(y)=ln(a)+b*ln(x)
	// A = ln(a) ; B = b

	vector <double> X = log(x); // X = ln(x)
	vector <double> Y = log(y); // Y = ln(y)

	vector <double> parameters = solveLSMApproximation(X, Y); // [a0, a1, a2, ..., am]
	double standar_deviation = calculateStandarDeviation(X, Y, parameters);

	cout << "\nStandar deviation = " << standar_deviation << "\n";

	parameters[0] = exp(parameters[0]); // a = exp(A)

	return parameters;
}
vector <double> solveLSMApproximation(vector <double> x, vector <double> y) {
	if (x.size() != y.size()) {
		throw invalid_argument("Vectors must have the same size.");
	}

	int N = x.size(); // Number of measurements (points)
	int m = 1; // Degree of the approximating polynomial

	vector <double> POWERX(2 * m, 0.0); // Sums of X coordinates in degrees from 1 to 2m
	for (int k = 1; k <= 2 * m; k++) {
		POWERX[k - 1] = sum(pow(x, k));
	}

	vector <vector<double>> SUMX(m + 1, vector<double>(m + 1, 0.0)); // Coefficient matrix, or matrix A in SLOUGH
	vector<double> PRAW(m + 1, 0.0); // The right part, or matrix B in SLOUGH
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= m; j++) {
			if (i != 0 || j != 0) {
				int k = i + j - 1;
				SUMX[i][j] = POWERX[k];
			}
			else {
				SUMX[i][j] = N;
			}
		}
		PRAW[i] = y * pow(x, i);
	}

	vector<double> parameters = solveGaussMethod(SUMX, PRAW); // [a0, a1, a2, ..., am]
	return parameters;
}
double calculateStandarDeviation(vector <double> x, vector <double> y, vector <double> parameters) {
	double standar_deviation = 0.0;

	int N = x.size(); // Number of measurements (points)
	int m = parameters.size() - 1; // Degree of the approximating polynomial

	for (int i = 0; i < N; i++) {
		double sum = 0.0;
		for (int j = 0; j <= m; j++) {
			sum += parameters[j] * pow(x[i], j); // (a0+a1x^1+a2x^2+...+amx^m) = Pm(x)
		}
		standar_deviation += pow(y[i] - sum, 2); // (y(x) - Pm(x))^2
	}

	standar_deviation /= (N - m - 1);
	return standar_deviation;
}