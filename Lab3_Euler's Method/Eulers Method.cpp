#include "Eulers Method.h"

vector<double> solveExplicitEulerMethod(vector<double> u, pair<double, double> t_boundaries, double omega, bool showSteps) {
	const double E = max(abs(u)) * 0.01; // Local precision, for stability conditions
	const double h_max = abs((t_boundaries.second - t_boundaries.first) * 0.001); // Max allowable integration step

	vector<double> y = u; // Solution vector
	vector<double> f(u.size());

	double t = t_boundaries.first; // Parametr t
	double h = 0.0; // Integration Step
	int iter = 0;

	while (t < t_boundaries.second) {
		iter++;
		f = calculateVectorF(y, t, omega);
		h = min(calculateVectorH(f, h_max, E));
		y = calculateVectorY(y, f, h);
		t += h;
		if (showSteps)
			cout << iter << "|\t" << "y| " << y << "\t\t" << "t| " << t << "\n";
	};
	return y;
}
vector<double> solveImplicitEulerMethod(vector<double> u, pair<double, double> t_boundaries, vector<double> lambda, bool showSteps) {
	const vector<double> E_acceptable = abs(u) * 0.01; // Precisions
	vector<double> E(3);

	const double h_min = abs((t_boundaries.second - t_boundaries.first) * 0.001); // Min allowable integration step
	const double h_max = abs((t_boundaries.second - t_boundaries.first) * 0.1); // Max allowable integration step
	const vector<vector<double>> A = fillMatrixA(lambda);
	const vector<double> b = fillVectorB(lambda);

	vector<double> y = u; // Solution vector
	vector<double> y_prev = u;
	vector<double> y_next = u;

	double h = h_min; // Integration Step
	double h_prev = h_min;
	double h_next = 0.0;

	double t = t_boundaries.first; // Parametr t
	double t_next = 0.0;
	int iter = 0;

	while (t < t_boundaries.second) {
		iter++;
		t_next = t + h;
		y_next = solveNewtonsMethod(y_next, y, h, A, b, precisions::E1, precisions::E2, 0.05);
		E = calculateLocalError(y_prev, y, y_next, h, h_prev);

		if (isLocalErrorExceedsPermissibleOne(E, E_acceptable)) {
			h = h / 2;
			t_next = t;
			y_next = y;
		}
		else {
			h_next = calculateNextStepOfIntegration(E, E_acceptable, h);
			if (h_next > h_max) {
				h_next = h_max;
			}
			if (showSteps)
				cout << iter << "|\t" << "y| " << y_next << "\t\t" << "t| " << t_next << "\n";
			y_prev = y;
			y = y_next;
			h_prev = h;
			h = h_next;
			t = t_next;
		}
	}
	return y;
}

double u0dt(vector<double> u, double t) { 
	if (t == 0)
		return 0;
	else
		return -1 * u[0] * u[1] + sin(t) / t; }
double u1dt(vector<double> u, double t, double x) { 
	return -1 * pow(u[1], 2) + (((2.5 + x / 40) * t) / (1 + pow(t, 2)));
}
vector<double> calculateVectorF(vector<double> y, double t, double x) {
	vector<double> f(y.size());
	f[0] = u0dt(y, t);
	f[1] = u1dt(y, t, x);
	return f; 
}
vector<double> calculateVectorH(vector<double> f, double h_max, double E) {
	vector<double> h(f.size());	
	for (int i = 0; i < h.size(); i++) {
		h[i] = E / (abs(f[i]) + (E / h_max));
	}
	return h;
}
vector<double> calculateVectorY(vector<double> y, vector<double> f, double h) {
	return y + f * h;
}

vector<vector<double>> fillMatrixA(vector<double> lambda) {
	vector<vector<double>> A(3, vector<double>(3));
	double l1 = lambda[0]; 
	double l2 = lambda[1]; 
	double l3 = lambda[2];
	A[0] = { 2 * l1 + 4 * l2 , 2 * (l1 - l2) , 2 * (l1 - l2) };
	A[1] = { 2 * (l1 - l2) , 2 * l1 + l2 + 3 * l3 , 2 * l1 + l2 - 3 * l3 };
	A[2] = { 2 * (l1 - l2) , 2 * l1 + l2 - 3 * l3 , 2 * l1 + l2 + 3 * l3 };
	A = A * (1.0 / 6);
	return A;
}
vector<double> fillVectorB(vector<double> lambda) {
	vector<double> b(3);
	double l1 = lambda[0]; 
	double l2 = lambda[1]; 
	double l3 = lambda[2];
	b = { 4 * l1 + 2 * l2 , 4 * l1 - l2 - 9 * l3 , 4 * l1 - l2 + 9 * l3 };
	b = b * (-1.0 / 6);
	return b;
}
vector<double> calculateLocalError(vector<double> y_prev, vector<double> y, vector<double> y_next,  double h, double h_prev) {
	vector<double> error(3);
	error = y_next - y - (y - y_prev) * (h / h_prev);
	error = error * (-h / (h + h_prev));
	return error; 
}
bool isLocalErrorExceedsPermissibleOne(vector<double> E, const vector<double> E_acceptable) {
	for (int i = 0; i < E.size();i++) {
		if (abs(E[i]) > E_acceptable[i]) {
			return true;
		}
	}
	return false; 
}
double calculateNextStepOfIntegration(vector<double> E, const vector<double> E_acceptable, double h) {
	return min(sqrt(E_acceptable / abs(E)) * h);
}