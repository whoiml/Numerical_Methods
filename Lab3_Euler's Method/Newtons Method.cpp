#include "Newtons Method.h"

vector<double> udt(vector<double> u, vector<vector<double>> A, vector<double> b) {
    return A * u - b;
}
vector<double> f(vector<double> y_next, vector<double> y, double h, vector<vector<double>> A, vector<double> b) {
    return (y_next - y) - (udt(y_next, A, b) * h);
}
vector<vector<double>> calculateMatrixJacobi(vector<double> x, vector<double> y, double h, const double M, const double E1, vector<vector<double>> A, vector<double> b) {
    vector<double> y_next = x;
    vector<double> ram_left(3);
    vector<double> ram_right(3);
    vector<vector<double>> Jacobi(3, vector<double>(3));

    double dy = 0.0;
    for (int i = 0; i < Jacobi.size(); i++) {
        dy = M * x[i];
        if (abs(dy) < E1) 
            dy = E1; 
        y_next[i] += dy;
        ram_left = f(y_next, y, h, A, b);
        ram_right = f(x, y, h, A, b);
        for (int j = 0; j < Jacobi.size(); j++) {
            Jacobi[j][i] = (ram_left[j] - ram_right[j]) / dy; 
        }
        y_next[i] -= dy;
    }
    return Jacobi; 
}

vector<double> solveNewtonsMethod(vector<double> x, vector<double> y, double h, vector<vector<double>> A, vector<double> b, const double E1, const double E2, const double M) {
    double delta1 = max(abs(f(x,y,h,A,b)));
    double delta2 = 1;

    vector<double> F(x.size());
    vector<double> solve_gauss(x.size());
    vector<double> ram_for_d2(x.size());
    vector<vector<double>> J(x.size(), vector<double>(x.size()));

    while ((delta1 > E1 || delta2 > E2)) {
        F = f(x, y, h, A, b);
        J = calculateMatrixJacobi(x, y, h, M, E1, A, b);
        solve_gauss = solveGaussMethod(J, -F);

        ram_for_d2 = x;
        x = x + solve_gauss;
        delta1 = max(abs(F));

        for (int i = 0; i < ram_for_d2.size(); i++) {
            ram_for_d2[i] = abs(ram_for_d2[i]) < 1 ? abs(solve_gauss[i]) : abs(solve_gauss[i] / ram_for_d2[i]);
        }
        delta2 = max(ram_for_d2);
    }
    return x; 
}