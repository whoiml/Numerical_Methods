#include "Eulers Method.h"

int main()
{
	cout.precision(10);

	cout << "Show iteration steps?\nEnter 1 or 0 : ";
	bool showSteps; cin >> showSteps;

	pair<double, double> t_boundaries(0,1);

	// Задание 1 - Явный метод Эйлера для 1й задачи
	double omega = 30; 
	vector<double> u = { 0, -0.412 };
	vector<double> solve_explicit = solveExplicitEulerMethod(u, t_boundaries, omega, showSteps);
	cout << "\nSolution is via the Explicit Euler method : " << solve_explicit << "\n";

	// Задание 2 - Неявный метод Эйлера для 4й задачи
	vector<double> u_new = { 10, 22, 9 };
	vector<double> lambda = { -100, -100, -100 };
	vector<double> solve_implicit = solveImplicitEulerMethod(u_new, t_boundaries, lambda, showSteps);
	cout << "\nSolution is via the Explicit Euler method : " << solve_implicit << "\n";

	return 0;
}
