#include <iostream>
#include "MatrixPrinter.h"

using namespace std;

class GaussReturnType {
public:
    double **matrixL;
    int *indexes;

    explicit GaussReturnType(int n) {
        matrixL = new double *[n];
        for (int i = 0; i < n; i++) {
            matrixL[i] = new double[n];
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrixL[i][j] = 0;
            }
        }

        indexes = new int[n];
        for (int i = 0; i < n; i++) {
            indexes[i] = i;
        }

    }
};

void gauss_recur(double **matrix, double **matrixL, int *indexes, int n, int p = 0) {
    if (p < n - 1) {
        for (int i = p + 1; i < n; i++) {
            // częściowy wybór elementu podstawowego
            if (matrix[indexes[p]][p] == 0) {
                double max = matrix[indexes[p + 1]][p];
                int max_idx = p + 1;
                for (int k = p + 1; k < n; k++) {
                    if (fabs(matrix[indexes[k]][p]) > max) {
                        max = fabs(matrix[indexes[k]][p]);
                        max_idx = k;
                    }
                }
                int tmp = indexes[p];
                indexes[p] = indexes[max_idx];
                indexes[max_idx] = tmp;
            }

            // wyznaczenie współczynnika
            double coeff = matrix[indexes[i]][p] / matrix[indexes[p]][p];
            matrixL[indexes[i]][p] = coeff;

            // eliminacja
            for (int j = p; j < n; j++) {
                matrix[indexes[i]][j] -= matrix[indexes[p]][j] * coeff;
            }
        }
        gauss_recur(matrix, matrixL, indexes, n, p + 1);
    }
}

GaussReturnType gauss(double **matrix, int n) {
    auto result = GaussReturnType(n);
    gauss_recur(matrix, result.matrixL, result.indexes, n);

    // wypełnienie macierzy L jednykami na przekątnej
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                result.matrixL[result.indexes[i]][j] = 1;
            }
        }
    }

    return result;
}

void vectorB(double *vector, double **matrixL, const int *indexes, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            vector[indexes[j]] -= vector[indexes[i]] * matrixL[indexes[j]][i];
        }
    }
}

double *solve(double **matrix, const double *vector, const int *indexes, int n) {
    auto *result = new double[n];
    double sum;
    for (int i = n - 1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += matrix[indexes[i]][j] * result[j];
        }
        result[i] = (vector[indexes[i]] - sum) / matrix[indexes[i]][i];
    }
    return result;
}

int main() {

    int n = 4;
    auto **A = new double *[n];
    A[0] = new double[n]{1, -20, 30, -4};
    A[1] = new double[n]{2, -40, -6, 50};
    A[2] = new double[n]{9, -180, 11, -12};
    A[3] = new double[n]{-16, 15, -140, 13};
    auto *b = new double[n]{35, 104, -366, -354};

    auto gaussReturnType = gauss(A, n);
    vectorB(b, gaussReturnType.matrixL, gaussReturnType.indexes, n);
    auto res = solve(A, b, gaussReturnType.indexes, n);

    MatrixPrinter::printMatrix(A, nullptr, n, "Macierz A");
    MatrixPrinter::printMatrix(A, gaussReturnType.indexes, n, "Macierz U");
    MatrixPrinter::printMatrix(gaussReturnType.matrixL, gaussReturnType.indexes, n, "Macierz L");
    MatrixPrinter::printVector(b, gaussReturnType.indexes, n, "Wektor b");
    MatrixPrinter::printVector(res, gaussReturnType.indexes, n, "Wektor X");

}