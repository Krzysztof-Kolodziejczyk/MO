#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

int MAX_ITER = 100;
double TOL_F = 1.e-5;
double TOL_X = 1.e-5;
int precission = 16;

void saveResult(double *x, int n, double residuum, double error, int iter, ofstream &file) {
    file << setprecision(precission) << iter << ")\t";
    for (int i = 0; i < n; i++) {
        file << setw(20) << x[i] << "  ";
    }
    file << "\t residuum = " << setw(20) << residuum << "       error = " << setw(20) << error << "\n";
}

void printVector(double *v, int n) {
    cout << setprecision(10);
    for (int i = 0; i < n; i++) {
        cout << v[i] << "  ";
    }
}

double calculateNorm(double *v, int n) {
    double max = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(v[i]) > max)
            max = fabs(v[i]);
    }
    return max;
}

double residuum(double **A, const double *x, const double *b, int n) {
    auto res = new double[n];
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        res[i] = sum - b[i];
    }
    return calculateNorm(res, n);
}

double errorEst(const double *v1, const double *v2, int n) {
    auto res = new double[n];
    for (int i = 0; i < n; i++)
        res[i] = v1[i] - v2[i];
    return calculateNorm(res, n);
}

double sum(double **A, const double *x, int i, int m) {
    double sum = 0.;
    for (int j = 0; j < m; j++) {
        if (i != j)
            sum += A[i][j] * x[j];
    }
    return sum;
}

double *Jacobi(double **A, const double *b, double *x, int n) {
    auto next_x = new double[n];
    int iter;
    double resi, error;
    ofstream file;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_7\jacobi)");
    for (iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++)
            next_x[i] = (b[i] - sum(A, x, i, n)) / A[i][i];
        resi = residuum(A, next_x, b, n);
        error = errorEst(x, next_x, n);
        saveResult(next_x, n, resi, error, iter, file);
        if (resi < TOL_F && error < TOL_X)
            break;
        for (int i = 0; i < n; i++)
            x[i] = next_x[i];
    }
    file.close();
    return x;
}


double *gSeidel(double **A, const double *b, double *x, int n) {
    auto *prev_x = new double[n];
    int iter;
    double resi, error;
    ofstream file;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_7\gSeidel)");
    for (iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++) {
            prev_x[i] = x[i];
            x[i] = (b[i] - sum(A, x, i, n)) / A[i][i];
        }
        resi = residuum(A, x, b, n);
        error = errorEst(x, prev_x, n);
        if (resi < TOL_F && error < TOL_X)
            break;
        saveResult(x, n, resi, error, iter, file);
    }
    file.close();
    return x;
}

double *sor(double **A, double *b, double *x, double omega, int n) {
    auto *prev_x = new double[n];
    int iter;
    double resi, error;
    ofstream file;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_7\sor)");
    for (iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++) {
            prev_x[i] = x[i];
            x[i] = (1. - omega) * x[i] + (omega / A[i][i]) * (b[i] - sum(A, x, i, n));
        }
        resi = residuum(A, x, b, n);
        error = errorEst(x, prev_x, n);
        if (resi < TOL_F && error < TOL_X)
            break;
        saveResult(x, n, resi, error, iter, file);
    }
    file.close();
    return x;
}

int main() {
    int n = 4;
    auto **A = new double *[n];
    A[0] = new double[n]{100., -1., 2., -3.};
    A[1] = new double[n]{1., 200., -4., 5.};
    A[2] = new double[n]{-2., 4., 300., -6.};
    A[3] = new double[n]{3., -5., 6., 400.};
    auto *b = new double[n]{116., -226., 912., -1174.};
    auto *init_x = new double[n]{2., 2., 2., 2.};
    printVector(Jacobi(A, b, init_x, n), n);

    init_x = new double[n]{2., 2., 2., 2.};
    cout << endl;
    printVector(gSeidel(A, b, init_x, n), n);

    init_x = new double[n]{2., 2., 2., 2.};
    cout << endl;
    printVector(sor(A, b, init_x, 0.5, n), n);
}