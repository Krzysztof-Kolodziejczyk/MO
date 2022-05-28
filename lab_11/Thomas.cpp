

#include "Thomas.h"

#include "iostream"

 using namespace std;

 double **Thomas::thomas(double **A, int n) {
    auto res = new double *[n];
    for (int i = 0; i < n; i++) {
        res[i] = new double[3];
    }
    res[0][0] = A[0][0];
    res[0][1] = A[0][1];
    res[0][2] = A[0][2];
    for (int i = 1; i < n; i++) {
        res[i][1] = A[i][1] - A[i][0] * A[i - 1][2] / res[i - 1][1];
        res[i][2] = A[i][2];
        res[i][0] = A[i][0] / res[i - 1][1];
    }
    return res;
}

 void Thomas::vectorB(double **matrixThomas, double *b, int n) {
    for (int i = 1; i < n; i++) {
        b[i] -= b[i - 1] * matrixThomas[i][0];
    }
}

 double *Thomas::solve(double **matrixThomas, const double *b, int n) {
    auto res = new double[n];
    res[n - 1] = b[n - 1] / matrixThomas[n - 1][1];
    for (int i = n - 2; i >= 0; i--) {
        res[i] = (b[i] - matrixThomas[i][2] * res[i + 1]) / matrixThomas[i][1];
    }
    return res;
}
