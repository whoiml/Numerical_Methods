#include "Gauss Method.h"

vector<double> operator / (const vector<double>& v1, const vector<double>& v2) {
    vector<double> result(v2.size());

    for (int i = 0; i < v2.size(); i++) {
        result[i] = v1[i] / v2[i];
    }
    return result;
}
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2) {
    vector<double> result(v2.size());

    for (int i = 0; i < result.size(); i++) {
        result[i] = 0;
        for (int j = 0; j < result.size(); j++) {
            result[i] += m1[i][j] * v2[j];
        }
    }
    return result;
}
vector<double> operator - (const vector<double>& v) {
    vector<double> result = v;
    for (int i = 0; i < v.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}
ostream& operator << (ostream& os, const vector<vector<double>>& A)
{
    int n = A.size();

    cout << endl;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            os << A[i][j] << " ";
        }
        cout << endl;
    }
    return os;
}
ostream& operator << (ostream& os, const vector<double>& A)
{
    int n = A.size();

    os << "(" << A[0];
    for (int i = 1; i < n; i++) {
        os << "; " << A[i];
    }
    os << ")" << endl << endl;
    return os;
}

bool isVectorWithSameCoordinates(vector<double> v)
{
    double ram = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (ram == v[i])
            continue;
        else
            return false;
    }
    return true;
}
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    vector<vector<double>> A_b(n, vector<double>(n + 1, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (j == n)
                A_b[i][j] = b[i];
            else
                A_b[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (isVectorWithSameCoordinates(A_b[i] / A_b[j]))
                return true;
        }
    }

    return false;
}
vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    if (isMatrixWithProportionalityRows(A, b)) {
        exit(errors::MATRIX_WITH_PROPORTIONAL_ROWS);
    }

    forwardGaussMethod(A, b);

    vector<double> x(n);
    backwardGaussMethod(A, b, x);
    return x;
}
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int k = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[k][i])) {
                k = j;
            }
        }

        swap(A[i], A[k]);
        swap(b[i], b[k]);

        double div = A[i][i];
        for (int j = i; j < n; j++) {
            A[i][j] /= div;
        }
        b[i] /= div;

        for (int j = i + 1; j < n; j++) {
            double mult = A[j][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= mult * A[i][k];
            }
            b[j] -= mult * b[i];
        }
    }
}
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    int n = A.size();

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
}