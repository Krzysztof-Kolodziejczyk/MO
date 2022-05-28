#include<iostream>
#include <iomanip>
#include <cfloat>
#include <fstream>
#include "cmath"
#include "calerf.h"
#include "Thomas.h"

using namespace std;

// parametry zadania
double t_max = 2.;
double r = 1.;
double a = 10.;
double D = 1.;

// parametry dla KMB
int x_size_PMCN, t_size_PMCN;
double KMB_lambda = 0.4;

// parametry dla PMCN
double PMCN_lambda = 1.;

void print_matrix(double **mat, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << setw(3) << mat[i][j] << " ";
        }
        cout << endl;
    }
}

// analityczne rozwiązanie funkcji w punkcie x,t
double analytical_solution(double x, double t) {
    return 1. - (r / x) * calerf::ERFC_L((x - r) / (2. * sqrt(D * t)));
}

// boundary condition for x = r + a
double right_boundary_cond(double t) {
    return 1. - (r / (r + a)) * calerf::ERFC_L(a / (2.0 * sqrt(D * t)));
}

double left_boundary_cond(double t) {
    return 0.;
}

double starting_cond(double x) {
    return 1.;
}

double **init_matrix(int t_size, int x_size) {
    auto **res = new double *[t_size];
    for (int i = 0; i < t_size; i++) {
        res[i] = new double[x_size];
    }
    return res;
}

void init_matrix_with_conditions(double **m, int t_size, int x_size, double dt, double h) {
    // inicjacja warunków brzegowych
    double t_iterator = dt;
    for (int i = 1; i < t_size; i++) {
        m[i][0] = left_boundary_cond(t_iterator);
        m[i][x_size - 1] = right_boundary_cond(t_iterator);
        t_iterator += dt;
    }

    // inicjajcja warunku początkowego
    double x_iterator = r;
    for (int i = 0; i < x_size; i++) {
        m[0][i] = starting_cond(x_iterator);
        x_iterator += h;
    }
}

double calculate_absolute_error(double **approx_matrix, double **analytical_matrix, int t_size, int x_size) {
    double max_error = 0.;
    double current_error;
    for (int i = 0; i < x_size; i++) {
        current_error = fabs(analytical_matrix[t_size - 1][i] - approx_matrix[t_size - 1][i]);
        if (current_error > max_error) {
            max_error = current_error;
        }
    }
    return max_error;
}

void save(ofstream &file, double log_h, double log_error) {
    file << setprecision(20) << log_h << "\t" << log_error << endl;
}

// sekcja KMB
double **analytical_solution_KMB_matrix(int t_size, int x_size, double dt, double h) {
    auto res = init_matrix(t_size, x_size);
    double t_iterator = 0.;
    double x_iterator;
    for (int it = 0; it < t_size_PMCN; it++) {
        x_iterator = r;
        for (int ih = 0; ih < x_size_PMCN; ih++) {
            res[it][ih] = analytical_solution(x_iterator, t_iterator);
            x_iterator += h;
        }
        t_iterator += dt;
    }
    return res;
}

void KMB_method(double **matrix, double h) {
    double x_iterator;
    for (int it = 1; it < t_size_PMCN; it++) {
        x_iterator = r + h;
        for (int ih = 1; ih < x_size_PMCN - 1; ih++) {
            matrix[it][ih] = KMB_lambda * matrix[it - 1][ih - 1] * (1. - (h / x_iterator))
                             + matrix[it - 1][ih] * (1. - (2. * KMB_lambda))
                             + KMB_lambda * matrix[it - 1][ih + 1] * (1. + (h / x_iterator));
            x_iterator += h;
        }
    }
}

// sekcja MPCN
double **init_triDiagonal_matrix(int x_size) {
    auto triDiagonal = new double *[x_size];
    for (int i = 0; i < x_size; i++) {
        triDiagonal[i] = new double[3];
    }
    return triDiagonal;
}

void fill_triDiagonal_matrix(double **A, int x_size, double h) {
    // lewy warunek brzegowy
    A[0][0] = NULL;
    A[0][1] = 1.;
    A[0][2] = 0.;

    // środkowa część macierzy
    double x_iteration = r + h;
    for (int i = 1; i < x_size - 1; i++) {
        A[i][0] = (PMCN_lambda / 2.) * (-1. + h / x_iteration);
        A[i][1] = 1. + PMCN_lambda;
        A[i][2] = (PMCN_lambda / 2.) * (-1. - h / 2.);
        x_iteration += h;
    }

    // prawy warunek brzegowy
    A[x_size - 1][0] = 0.;
    A[x_size - 1][1] = 1.;
    A[x_size - 1][2] = NULL;
}

// k - poziom niewiadomych, k-1 poziom znanych wartosci
double *init_and_fill_B_vector(double **matrix, int x_size, double h, double t, int k) {
    auto *b = new double[x_size];

    // lewy warunek brzegowy
    b[0] = left_boundary_cond(t);

    // srodkowa częśc wektora
    double x_iteration = r + h;
    for (int i = 1; i < x_size - 1; i++) {
        b[i] = (PMCN_lambda / 2.) * matrix[k - 1][i - 1] * (1. - h / x_iteration)
               + matrix[k-1][i] * (1. - PMCN_lambda)
               + (PMCN_lambda / 2.) * matrix[k - 1][i + 1] * (1. + h / x_iteration);
        x_iteration += h;
    }

    // prawy warunek brzegowy
    b[x_size - 1] = right_boundary_cond(t);

    return b;
}

void PMCN_method(double **matrix, double h, double dt, int t_size, int x_size) {
    auto A = init_triDiagonal_matrix(x_size);
    fill_triDiagonal_matrix(A, x_size, h);
    Thomas::thomas(A, x_size);

    double *b;
    double t_iteration = dt;
    for (int k = 1; k < t_size; k++) {
        b = init_and_fill_B_vector(matrix, x_size, h, t_iteration, k);
        Thomas::vectorB(A, b, x_size);
        auto u_result = Thomas::solve(A, b, x_size);
        for (int i = 1; i < x_size - 1; i++) {
            matrix[k][i] = u_result[i];
        }
        t_iteration += dt;
    }

}


// MAIN METHODS
void KMB(){
    double dt;
    ofstream KMB_error_file;
    KMB_error_file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\KMB\KMB_error.txt)");
    double **KMB_matrix;
    double **KMB_analytical_matrix;

    for (double h = 0.5; h > 0.01; h -= 0.01) {
        dt = (h * h) * KMB_lambda;

        t_size_PMCN = floor(t_max / dt);
        x_size_PMCN = floor(a / h);

        KMB_matrix = init_matrix(t_size_PMCN, x_size_PMCN);
        init_matrix_with_conditions(KMB_matrix, t_size_PMCN, x_size_PMCN, dt, h);

        KMB_method(KMB_matrix, h);

        KMB_analytical_matrix = analytical_solution_KMB_matrix(t_size_PMCN, x_size_PMCN, dt, h);

        auto kmb_error = calculate_absolute_error(KMB_matrix, KMB_analytical_matrix, t_size_PMCN, x_size_PMCN);

        save(KMB_error_file, log10(h), log10(kmb_error));

    }

    KMB_error_file.close();
}

void PMCN(){
    double dt;
    ofstream PMCN_error_file;
    PMCN_error_file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN\PMCN_error.txt)");
    double **PMCN_matrix;
    double **PMCN_analytical_matrix;

    for (double h = 0.5; h > 0.01; h -= 0.01) {
        dt = (h * h) * PMCN_lambda;

        t_size_PMCN = floor(t_max / dt);
        x_size_PMCN = floor(a / h);

        PMCN_matrix = init_matrix(t_size_PMCN, x_size_PMCN);
        init_matrix_with_conditions(PMCN_matrix, t_size_PMCN, x_size_PMCN, dt, h);

        PMCN_method(PMCN_matrix, h,dt, t_size_PMCN, x_size_PMCN);

        PMCN_analytical_matrix = analytical_solution_KMB_matrix(t_size_PMCN, x_size_PMCN, dt, h);

        auto kmb_error = calculate_absolute_error(PMCN_matrix, PMCN_analytical_matrix, t_size_PMCN, x_size_PMCN);

        save(PMCN_error_file, log10(h), log10(kmb_error));

    }

    PMCN_error_file.close();
}


int main() {
    PMCN();
}

