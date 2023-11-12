#include "Methods.h"

int main()
{
	double a = 0.0, b = PI;
	double A = -1, B = 1;
	double C = -1, D = 1;

	cout << "\n\nE = " << E4 << "\n______________________";
	cout << "\nTrapezoid Method\n"; solveTrapezoidMethod(a, b, f_v19, E4);
	cout << "\nSimpson's method\n"; solveSimpsonsMethod(a, b, f_v19, E4);
	cout << "\nSimpson's Cubature Method\n"; solveSimpsonsCubatureMethod(A, B, C, D, f_v32, E4);
	
	cout << "\n\nE = " << E5 << "\n______________________";
	cout << "\nTrapezoid Method\n"; solveTrapezoidMethod(a, b, f_v19, E5);
	cout << "\nSimpson's method\n"; solveSimpsonsMethod(a, b, f_v19, E5);
	cout << "\nSimpson's Cubature Method\n"; solveSimpsonsCubatureMethod(A, B, C, D, f_v32, E5);
	return 0;
}