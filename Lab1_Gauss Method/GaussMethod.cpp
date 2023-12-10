#include "GaussMethod.h"

template<typename T>
vector<T> operator / (vector<T>& v1, vector<T>& v2) {
    vector<T> result(v2.size());

    for (int i = 0; i < v2.size(); i++) {
        result[i] = v1[i] / v2[i];
    }
    return result;
}

template<typename T>
vector<T> operator * (vector<vector<T>>& m1, vector<T>& v2) {
    vector<T> result(v2.size());

    for (int i = 0; i < result.size(); i++) {
        result[i] = 0;
        for (int j = 0; j < result.size(); j++) {
            result[i] += m1[i][j] * v2[j];
        }
    }
    return result;
}

template<typename T>
ostream& operator << (ostream& os, vector<vector<T>>& A) {
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

template<typename T>
ostream& operator << (ostream& os, vector<T>& A) {
    int n = A.size();

    os << "(" << A[0];
    for (int i = 1; i < n; i++) {
        os << "; " << A[i];
    }
    os << ")" << endl << endl;
    return os;
}

template<typename T>
istream& operator >> (istream& is, vector<vector<T>>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            is >> A[i][j];
        }
    }
    return is;
}

template<typename T>
istream& operator >> (istream& is, vector<T>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        is >> A[i];
    }
    return is;
}

template<typename T>
bool isVectorWithSameCoordinates(vector<T> v) {
    T ram = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (ram == v[i])
            continue;
        else
            return false;
    }
    return true;
}

template<typename T>
bool isMatrixWithProportionalityRows(vector<vector<T>> A, vector<T> b) {
    int n = A.size();

    vector<vector<T>> A_b(n, vector<T>(n + 1, 0));
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

template<typename T>
vector<T> solveGaussMethod(vector<vector<T>> A, vector<T> b) {
    int n = A.size();

    if (isMatrixWithProportionalityRows(A, b)) {
        exit(errors::MATRIX_WITH_PROPORTIONAL_ROWS);
    }

    forwardGaussMethod(A, b);

    vector<T> x(n);
    backwardGaussMethod(A, b, x);
    return x; 
}

template<typename T>
void forwardGaussMethod(vector<vector<T>>& A, vector<T>& b) {
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

        T div = A[i][i];
        for (int j = i; j < n; j++) {
            A[i][j] /= div;
        }
        b[i] /= div;

        for (int j = i + 1; j < n; j++) {
            T mult = A[j][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= mult * A[i][k];
            }
            b[j] -= mult * b[i];
        }
    }
}

template<typename T>
void backwardGaussMethod(vector<vector<T>> A, vector<T> b, vector<T>& x) {
    int n = A.size();

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
}

template<typename T>
vector<T> calculateResidualVector(vector<vector<T>> A, vector<T> b, vector<T> x) {
    int n = A.size();
    double Ax = 0;

    vector<T> F(n);
    for (int i = 0; i < n; i++) {
        Ax = 0;
        for (int j = 0; j < n; j++) {
            Ax += A[i][j] * x[j];
        }
        F[i] = Ax - b[i];
    }
    return F;
}

template<typename T>
T calculateNormaOfRV(vector<T> F) {
    int n = F.size();
    T norma = 0;

    for (int i = 0; i < n; i++) {
        norma += pow(F[i], 2);
    }
    norma = sqrt(norma);
    return norma;
}

template<typename T>
T calculateCalcError(vector<T> x1, vector<T> x2) {
    int n = x1.size();
    T calc_error = 0;

    T max1 = 0, max2 = 0;
    for (int i = 0; i < n; i++) {
        if (x2[i] - x1[i] > max1)
            max1 = x2[i] - x1[i];
        if (x1[i] > max2)
            max2 = x1[i];
    }

    calc_error = max1 / max2;
    return calc_error;
}

template<typename T>
vector<T> solveFactorizationLDLT(vector<vector<T>> A, vector<T> b) {
    int n = A.size();

    vector<T> D(n);
    vector<vector<T>> L(n, vector<T>(n, 0.0));
    vector<vector<T>> LT(n, vector<T>(n, 0.0));

    // Units on the main diagonal
	for (int i = 0; i < 3; ++i) {
		L[i][i] = 1;
		LT[i][i] = 1;
	}

    // Iteration cycle by columns
    for (int j = 0; j < n; j++) {
        // The amount we take away
        T sum = 0.0;
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
    vector<T> y(n);
    vector<T> z(n);
    for (int i = 0; i < n; i++) {
        T sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;

        z[i] = y[i] / D[i];
    }

	// LTx = z
    vector<T> x(n);
    for (int i = n - 1; i >= 0; i--) {
        T sum = 0.0;
        for (int j = n - 1; j > i; j--) {
            sum += LT[i][j] * x[j];
        }
        x[i] = z[i] - sum;
    }
    return x;
}
