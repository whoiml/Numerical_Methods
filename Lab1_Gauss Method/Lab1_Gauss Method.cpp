#include "GaussMethod.h"

int main()
{
    cout << "Choose how you want to use the Gauss method :\n"
            "1 - Enter the matrices yourself\n"
            "2 - Use the matrix according to the task condition\n"
            "Choise : ";
    int choise = 0; cin >> choise;
    while (choise < 1 || choise > 2) {
        cout << "Choose one of the options.";
        cin >> choise;
    }

    vector<vector<double>> A;
    vector<double> b;
    int N = 0;

    switch (choise) {
    case 1:
        cout << "\nEnter the size of the matrix A : ";
        cin >> N;
        while (N < 1) {
            cout << "Enter the correct size of the matrix.";
            cin >> N;
        }
        A.assign(N,vector<double>(N));
        b.resize(N);
        cout << "Enter matrix A :\n"; cin >> A;
        cout << "Enter matrix b :\n"; cin >> b;
        break; 
    case 2 :
        A = { { 1, 1, 1, 1 },
            { 1, 2, -2, 3 },
            { 2, 0, 1, 0 },
            { 3, 1, 2, 2 } };
        b = { 10, 11, 5, 19 };
        break;
    }

    // Задание 1 - Решение СЛАУ методом Гауса
    vector<double> x1 = solveGaussMethod(A, b);

    cout << "SLOUGH's Solution 1: " << endl << "x1 = " << x1;

    // Задание 2 - Вектор невязки и его норма
    vector<double> F = calculateResidualVector(A, b, x1); 
    double norma_of_RV = calculateNormaOfRV(F);

    cout << "The discrepancy vector: " << endl << "F = " << F;
    cout << "Norma of RV = " << norma_of_RV << endl << endl;

    // Задание 3 - Поиск относительной погрешности
    vector<double> Ax = A * x1;
    vector<double> x2 = solveGaussMethod(A, Ax);
    double calc_error = calculateCalcError(x1, x2);

    cout << "SLOUGH's Solution 2: " << endl << "x2 = " << x2;
    cout << "Calculation error = " << calc_error << endl << endl;

    // Задание 4 - LDLT факторизация
    double l1 = 1, l2 = 1000, l3 = 1000000;    
    vector<vector<double>> A_LDLT = { {2 * l1 + 4 * l2, 2 * (l1 - l2), 2 * (l1 - l2)},
                                 {2 * (l1 - l2), 2 * l1 + l2 + 3 * l3, 2 * l1 + l2 - 3 * l3},
                                 {2 * (l1 - l2), 2 * l1 + l2 - 3 * l3, 2 * l1 + l2 + 3 * l3 } };
    vector<double> b_LDLT = { -4 * l1 - 2 * l2, -4 * l1 + l2 + 9 * l3, -4 * l1 + l2 - 9 * l3 };

    vector<double> x3 = solveFactorizationLDLT(A_LDLT, b_LDLT);
    cout << "LDLT Solution : " << endl << "x3 = " << x3;
    return 0;
}

// Matrix by condition
//vector<vector<double>> A = { { 1, 1, 1, 1 },
//                             { 1, 2, -2, 3 },
//                             { 2, 0, 1, 0 },
//                             { 3, 1, 2, 2 } };
//vector<double> b = { 10, 11, 5, 19 };

// Matrix with proportionality rows
//vector<vector<double>> A = { { 1, 1, 1, 1 },
//                 { 2, 2, 2, 2 },
//                 { 1, 2, 3, 3 },
//                 { 1, 2, 3, 4 } };
//vector<double> b = { 10, 20, 19, 30 };
