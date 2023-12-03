#include "Vector.h"

vector<double> operator + (const vector<double>& v1, const vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw invalid_argument("Vectors must have the same size.");
    }

    vector<double> result(v1.size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}
vector<double> operator - (const vector<double>& v1, const vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw invalid_argument("Vectors must have the same size.");
    }

    vector<double> result(v1.size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = v1[i] - v2[i];
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
vector<double> operator * (const vector<double>& v, double scalar) {
    vector<double> result(v.size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = v[i] * scalar;
    }
    return result;
}
vector<double> operator / (const vector<double>& v1, const vector<double>& v2) {
    vector<double> result(v2.size());

    for (int i = 0; i < v2.size(); i++) {
        result[i] = v1[i] / v2[i];
    }
    return result;
}
vector<vector<double>> operator * (const vector<vector<double>>& m, double scalar) {
    vector<vector<double>> result(m.size(), vector<double>(m.size()));
    for (int i = 0; i < result.size(); i++) {
        result[i] = m[i] * scalar;
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
    os << ")";
    return os;
}

vector<double> abs(vector<double> v) {
    vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++) {
        result[i] = abs(v[i]);
    }
    return result;
}
vector<double> sqrt(vector<double> v) {
    vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++) {
        result[i] = sqrt(v[i]);
    }
    return result;
}
double max(vector<double> v) {
    double max = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (v[i] > max) {
            max = v[i];
        }
    }
    return max;
}
double min(vector<double> v) {
    double min = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (v[i] < min) {
            min = v[i];
        }
    }
    return min;
}