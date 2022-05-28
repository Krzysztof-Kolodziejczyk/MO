

#include "LU.h"

#include <cmath>

 void LU::gauss_recur(double **A, int *indexes, int n, int p) {
    if (p < n - 1) {
        for (int i = p + 1; i < n; i++) {
            // częściowy wybór elementu podstawowego
            if (A[indexes[p]][p] == 0) {
                double max = A[indexes[p + 1]][p];
                int max_idx = p + 1;
                for (int k = p + 1; k < n; k++) {
                    if (fabs(A[indexes[k]][p]) > max) {
                        max = fabs(A[indexes[k]][p]);
                        max_idx = k;
                    }
                }
                int tmp = indexes[p];
                indexes[p] = indexes[max_idx];
                indexes[max_idx] = tmp;
            }

            // wyznaczenie współczynnika
            double coeff = A[indexes[i]][p] / A[indexes[p]][p];
            // eliminacja
            for (int j = p; j < n; j++) {
                A[indexes[i]][j] -= A[indexes[p]][j] * coeff;
            }
            A[indexes[i]][p] = coeff;
        }
        gauss_recur(A, indexes, n, p + 1);
    }
}

int *LU::gauss(double **A, int n) {
    int *indexes = new int[n];
    for (int i = 0; i < n; i++) {
        indexes[i] = i;
    }
    LU::gauss_recur(A, indexes, n);

    return indexes;
}

void LU::vectorB(double *b, double **A, const int *indexes, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            b[indexes[j]] -= b[indexes[i]] * A[indexes[j]][i];
        }
    }
}

double *LU::solve(double **A, const double *b, const int *indexes, int n) {
    auto *result = new double[n];
    double sum;
    for (int i = n - 1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[indexes[i]][j] * result[j];
        }
        result[i] = (b[indexes[i]] - sum) / A[indexes[i]][i];
    }
    return result;
}
