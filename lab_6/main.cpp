#include "MatrixPrinter.h"

using namespace std;

double **thomas(double **A, int n) {
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

void vectorB(double **matrixThomas, double *b, int n) {
    for (int i = 1; i < n; i++) {
        b[i] -= b[i - 1] * matrixThomas[i][0];
    }
}

double *solve(double **matrixThomas, const double *b, int n) {
    auto res = new double[n];
    res[n - 1] = b[n - 1] / matrixThomas[n - 1][1];
    for (int i = n - 2; i >= 0; i--) {
        res[i] = (b[i] - matrixThomas[i][2] * res[i + 1]) / matrixThomas[i][1];
    }
    return res;
}

int main() {
    int n = 6;
    auto **A = new double *[n];
    A[0] = new double[n]{NULL, 10.0, 1.0 / 2};
    A[1] = new double[n]{1.0 / 3.0, 20.0, 1.0 / 4.0};
    A[2] = new double[n]{1.0 / 5.0, 30.0, 1.0 / 6.0};
    A[3] = new double[n]{1.0 / 7.0, 30.0, 1.0 / 8.0};
    A[4] = new double[n]{1.0 / 9.0, 20.0, 1.0 / 10.0};
    A[5] = new double[n]{1.0 / 11.0, 10.0, NULL};
    auto *b = new double[n]{31.0, 165.0 / 4.0, 917.0 / 30.0, 851.0 / 28.0, 3637.0 / 90.0, 332.0 / 11.0};

    auto resultThomas = thomas(A, n);

    MatrixPrinter::printMatrix(A, nullptr, n, 3, "Macierz A");
    MatrixPrinter::printVector(b, nullptr, n, "Wektor B");
    MatrixPrinter::printMatrix(resultThomas, nullptr, n, 3, "Macierz Thomas");

    vectorB(resultThomas, b, n);
    MatrixPrinter::printVector(b, nullptr, n, "Wektor B po transformacji");

    auto vectorX = solve(resultThomas, b, n);
    MatrixPrinter::printVector(vectorX, nullptr, n, "Wektor X");
}