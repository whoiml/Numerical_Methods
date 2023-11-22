#include "Newtons Method.h"
#include "Gauss Method.h"

double f1(double x1, double x2) { return (sin(x1 + x2) - 1.3 * x1 - 0.1); }
double f1dx1(double x1, double x2) { return (cos(x1 + x2) - 1.3); }
double f1dx2(double x1, double x2) { return cos(x1 + x2); }
double f2(double x1, double x2) { return (pow(x1, 2) + pow(x2, 2) - 1); }
double f2dx1(double x1, double x2) { return 2 * x1; }
double f2dx2(double x1, double x2) { return 2 * x2; }

vector<vector<double>> calculateMatrixJ(double x1, double x2, double k) {
    return vector<vector<double>> { {(f1(x1 + x1 * k, x2) - f1(x1, x2)) / k / x1, (f1(x1, x2 + x2 * k) - f1(x1, x2)) / k / x2},
        { (f2(x1 + x1 * k, x2) - f2(x1, x2)) / k / x1, (f2(x1, x2 + x2 * k) - f2(x1, x2)) / k / x2 } };
}

void solveNewtonsMethod(double x1, double x2, const double E1, const double E2, const int max_iter, const double M) {
    double delta1 = max(abs(f1(x1, x2)), abs(f2(x1, x2)));
    double delta2 = 1;

    cout << "Initial approximation : ( " << x1 << " ; " << x2 << " )\n";

    /* M - Parameter is a relative increment.
    If it is NULL, then we consider the Jacobi matrix to be an analytical method,
    otherwise a numerically finite difference method, where it is needed. */
    if (M != NULL) { cout << "Relative increment = " << M << "\n"; }

    int iter = 0;
    while ((delta1 > E1 || delta2 > E2) && iter < max_iter) {
        iter++;

        cout << iter << "|\t" << "delta1| " << delta1 << "\t\t" << "delta2| " << delta2 << endl;

        vector<double> F{ f1(x1, x2), f2(x1, x2) };
        vector<vector<double>> J;

        if (M == NULL) {
            J = { {f1dx1(x1, x2), f1dx2(x1, x2)}, {f2dx1(x1, x2), f2dx2(x1, x2)} };
        }
        else {
            J = calculateMatrixJ(x1, x2, M);
        }

        vector<double> solution = solveGaussMethod(J, -F);
        x1 += solution[0];
        x2 += solution[1];

        delta1 = abs(F[0]);
        for (int i = 1; i < F.size(); i++) {
            if (delta1 < abs(F[i])) {
                delta1 = abs(F[i]);
            }
        }

        double v1 = abs(x1) < 1 ? abs(solution[0]) : abs(solution[0] / x1);
        double v2 = abs(x2) < 1 ? abs(solution[1]) : abs(solution[1] / x2);
        delta2 = max(v1, v2);
    }
    cout << "\n" << iter << "|\t" << "x1| " << x1 << "\t" << "x2| " << x2 << "\n";
    cout << "____________________________________________________________________________\n\n";
}
