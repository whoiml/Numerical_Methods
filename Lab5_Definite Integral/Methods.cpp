#include "Methods.h"

double f_v19(double x) { return 1 / (3 + 2 * cos(x)); }
double f_v32(int i, int j, int n)
{
	double A = -1, B = 1;
	double C = -1, D = 1;
	double h_x = (B - A) / (2 * n);
	double h_y = (D - C) / (2 * n);
	double x = A + i * h_x, y = C + j * h_y;
	if (x > 1 || y > 1)
		return 0;
	else
		return 4 - pow(x, 2) - pow(y, 2);
}

double solveTrapezoidMethod(double a, double b, double(*f)(double), const double E) {
	int n = 1;
	double h = b - a;
	double integral_prev = 0.0;
	double integral = h * (f(a) + f(b)) / 2.0;

	do {
		integral_prev = integral;
		double sum = 0.0;

		for (int i = 1; i <= n; ++i) {
			double x = a + i * h;
			sum += f(x);
		}

		integral = (integral + h * sum) / 2.0;
		h /= 2.0;
		n *= 2.0;
	} while (abs(integral - integral_prev) > E);

	double calc_error = abs((integral_prev - integral) / (pow(0.5, 2) - 1));

	cout << "- Result: " << integral << "\n";
	cout << "- Calculation error: " << calc_error << "\n";
	return integral;
}
double solveSimpsonsMethod(double a, double b, double(*f)(double), const double E) {
	int n = 1;
	double h = (b - a) / n;
	double integral_prev = 0;
	double integral = h * (f(a) + f(b)) / 2.0;

	do {
		integral_prev = integral;
		double sum = 0;

		for (int i = 1; i <= n; i++) {
			double x = a + (i - 0.5) * h;
			sum += f(x);
		}

		integral = (integral_prev + h * sum) / 2;
		h /= 2.0;
		n *= 2.0;
	} while (abs(integral - integral_prev) > 15 * E);

	double calc_error = abs((integral_prev - integral) / (pow(0.5, 4) - 1));

	cout << "- Result: " << integral << "\n";
	cout << "- Calculation error: " << calc_error << "\n";
	return integral;
}
double solveSimpsonsCubatureMethod(double a, double A, double b, double B, double (*f)(int, int, int), const double e)
{
	int n = 4;
	int m = 4;
	double h_x = (A - a) / (2 * n);
	double h_y = (B - b) / (2 * m);

	double integral = 0.0;
	double integral_prev = 0.0;

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			integral += f(2 * i, 2 * j, n) + f(2 * i + 2, 2 * j, n) + f(2 * i + 2, 2 * j + 2, n) + f(2 * i, 2 * j + 2, n) +
				4 * (f(2 * i + 1, 2 * j, n) + f(2 * i + 2, 2 * j + 1, n) + f(2 * i + 1, 2 * j + 2, n) + f(2 * i, 2 * j + 1, n)) +
				16 * f(2 * i + 1, 2 * j + 1, n);
		}
	}
	integral *= h_x * h_y / 9;

	do
	{
		n *= 2;
		h_x = (A - a) / (2 * n);
		h_y = (B - b) / (2 * m);
		integral_prev = integral;
		integral = 0;
		for (int i = 0; i <= n; i++)
		{
			for (int j = 0; j <= m; j++)
			{
				integral += f(2 * i, 2 * j, n) + f(2 * i + 2, 2 * j, n) + f(2 * i + 2, 2 * j + 2, n) + f(2 * i, 2 * j + 2, n) +
					4 * (f(2 * i + 1, 2 * j, n) + f(2 * i + 2, 2 * j + 1, n) + f(2 * i + 1, 2 * j + 2, n) + f(2 * i, 2 * j + 1, n)) +
					16 * f(2 * i + 1, 2 * j + 1, n);
			}
		}
		integral *= h_x * h_y / 9;
	} while (abs(integral - integral_prev) > 15 * e);

	cout << "- Result: " << integral << "\n";
	return integral;
}