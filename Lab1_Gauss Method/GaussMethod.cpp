#include "GaussMethod.h"

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

ostream& operator << (ostream& os, const vector<vector<double>>& A) {
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

ostream& operator << (ostream& os, const vector<double>& A) {
    int n = A.size();

    os << "(" << A[0];
    for (int i = 1; i < n; i++) {
        os << "; " << A[i];
    }
    os << ")" << endl << endl;
    return os;
}

istream& operator >> (istream& is, vector<vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            is >> A[i][j];
        }
    }
    return is;
}

istream& operator >> (istream& is, vector<double>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        is >> A[i];
    }
    return is;
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

vector<double> calculateResidualVector(vector<vector<double>> A, vector<double> b, vector<double> x) {
    int n = A.size();
    double Ax = 0;

    vector<double> F(n);
    for (int i = 0; i < n; i++) {
        Ax = 0;
        for (int j = 0; j < n; j++) {
            Ax += A[i][j] * x[j];
        }
        F[i] = Ax - b[i];
    }
    return F;
}

double calculateNormaOfRV(vector<double> F) {
    int n = F.size();
    double norma = 0;

    for (int i = 0; i < n; i++) {
        norma += pow(F[i], 2);
    }
    norma = sqrt(norma);
    return norma;
}

double calculateCalcError(vector<double> x1, vector<double> x2) {
    int n = x1.size();
    double calc_error = 0;

    double max1 = 0, max2 = 0;
    for (int i = 0; i < n; i++) {
        if (x2[i] - x1[i] > max1)
            max1 = x2[i] - x1[i];
        if (x1[i] > max2)
            max2 = x1[i];
    }

    calc_error = max1 / max2;
    return calc_error;
}

vector<double> solveLDLT_Factorization(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    vector<double> D(n);
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> LT(n, vector<double>(n, 0.0));

    // Units on the main diagonal
	for (int i = 0; i < 3; ++i) {
		L[i][i] = 1;
		LT[i][i] = 1;
	}

    // Iteration cycle by columns
    for (int j = 0; j < n; j++) {
        // The amount we take away
        double sum = 0.0;
        for (int i = 0; i < j; i++) {
            sum += D[i] * pow(L[j][i], 2);
        }

        // Fill matrix D
        D[j] = A[j][j] - sum;

        // Fill matrix L
        for (int i = j + 1; i < n; i++) {
            L[i][j] = (A[i][j] - sum) / D[j];
        }
    }

    // Fill matrix LT
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            LT[i][j] = L[j][i];
        }
    }

	// Ly = B ; Dz = y
    vector<double> y(n);
    vector<double> z(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;

        z[i] = y[i] / D[i];
    }

	// LTx = z
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = n - 1; j > i; j--) {
            sum += LT[i][j] * x[j];
        }
        x[i] = z[i] - sum;
    }
    return x;
}