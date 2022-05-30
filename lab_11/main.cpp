#include<iostream>
#include <iomanip>
#include <cfloat>
#include <fstream>
#include "cmath"
#include "calerf.h"
#include "Thomas.h"
#include "LU.h"

using namespace std;

// parametry zadania
double t_max = 2.;
double r = 1.;
double a = 10.;
double D = 1.;

// parametry dla KMB
int x_size, t_size;
double KMB_lambda = 0.4;

// parametry dla PMCN_THOMAS
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
    for (int i = 0; i < t_size; i++) {
        for (int j = 0; j < x_size; j++) {
            res[i][j] = 0.;
        }
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

double calculate_absolute_error(double **approx_matrix, double **analytical_matrix, int t_size, int x_size, int lvl) {
    double max_error = 0.;
    double current_error;
    for (int i = 0; i < x_size; i++) {
        current_error = fabs(analytical_matrix[i][lvl] - approx_matrix[i][lvl]);
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
    for (int it = 0; it < t_size; it++) {
        x_iterator = r;
        for (int ih = 0; ih < x_size; ih++) {
            res[it][ih] = analytical_solution(x_iterator, t_iterator);
            x_iterator += h;
        }
        t_iterator += dt;
    }
    return res;
}

void KMB_method(double **matrix, double h) {
    double x_iterator;
    for (int it = 1; it < t_size; it++) {
        x_iterator = r + h;
        for (int ih = 1; ih < x_size - 1; ih++) {
            matrix[it][ih] = KMB_lambda * matrix[it - 1][ih - 1] * (1. - (h / x_iterator))
                             + matrix[it - 1][ih] * (1. - (2. * KMB_lambda))
                             + KMB_lambda * matrix[it - 1][ih + 1] * (1. + (h / x_iterator));
            x_iterator += h;
        }
    }
}

// sekcja MPCN THOMAS
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
        A[i][0] = (PMCN_lambda / 2.) * (-1. + (h / x_iteration));
        A[i][1] = 1. + PMCN_lambda;
        A[i][2] = (PMCN_lambda / 2.) * (-1. - (h / x_iteration));
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
    b[0] = matrix[k - 1][0];

    // srodkowa częśc wektora
    double x_iteration = r + h;
    for (int i = 1; i < x_size - 1; i++) {
        b[i] = (PMCN_lambda / 2.) * matrix[k - 1][i - 1] * (1. - (h / x_iteration))
               + matrix[k - 1][i] * (1. - PMCN_lambda)
               + (PMCN_lambda / 2.) * matrix[k - 1][i + 1] * (1. + (h / x_iteration));
        x_iteration += h;
    }

    // prawy warunek brzegowy
    b[x_size - 1] = matrix[k - 1][x_size - 1];

    return b;
}

void PMCN_method(double **matrix, double h, double dt, int t_size, int x_size) {
    auto A = init_triDiagonal_matrix(x_size);
    fill_triDiagonal_matrix(A, x_size, h);
    auto matrixThomas = Thomas::thomas(A, x_size);

    double *b;
    double t_iteration = dt;
    for (int k = 1; k < t_size; k++) {
        b = init_and_fill_B_vector(matrix, x_size, h, t_iteration, k);
        Thomas::vectorB(matrixThomas, b, x_size);
        auto u_result = Thomas::solve(matrixThomas, b, x_size);
        for (int i = 1; i < x_size - 1; i++) {
            matrix[k][i] = u_result[i];
        }
        t_iteration += dt;
    }

}

// sekcja PMCN LU
double **init_matrix_for_LU(int x_size, double h) {
    double **res = init_matrix(x_size, x_size);

    // lewy warunek brzegowy
    res[0][0] = 1.;
    res[0][1] = 0.;

    // środkowa część macierzy
    double x_iteration = r + h;
    for (int i = 1; i < x_size - 1; i++) {
        res[i][i - 1] = (PMCN_lambda / 2.) * (-1. + (h / x_iteration));
        res[i][i] = 1. + PMCN_lambda;
        res[i][i + 1] = (PMCN_lambda / 2.) * (-1. - (h / x_iteration));
        x_iteration += h;
    }

    // prawy warunek brzegowy
    res[x_size - 1][x_size - 2] = 0.;
    res[x_size - 1][x_size - 1] = 1.;

    return res;
}

void PMCN_LU_method(double **matrix, double h, double dt, int t_size, int x_size) {
    auto A = init_matrix_for_LU(x_size, h);
    auto indexes = LU::gauss(A, x_size);

    double *b;
    double t_iteration = dt;
    for (int k = 1; k < t_size; k++) {
        b = init_and_fill_B_vector(matrix, x_size, h, t_iteration, k);
        LU::vectorB(b, A, indexes, x_size);
        auto u_result = LU::solve(A, b, indexes, x_size);
        for (int i = 1; i < x_size - 1; i++) {
            matrix[k][i] = u_result[i];
        }
        t_iteration += dt;
    }
}


// MAIN METHODS
void KMB(double start_h, double stop_h, double step, const string &file_path, bool error_case, bool analytical_case,
         bool time_case) {
    double dt;
    ofstream file;
    file.open(file_path);
    double **KMB_matrix;
    double **KMB_analytical_matrix;

    for (double h = start_h; h >= stop_h; h -= step) {

        // couting section
        dt = (h * h) * KMB_lambda;
        t_size = floor(t_max / dt);
        x_size = floor(a / h);
        KMB_matrix = init_matrix(t_size, x_size);
        init_matrix_with_conditions(KMB_matrix, t_size, x_size, dt, h);
        KMB_method(KMB_matrix, h);
        KMB_analytical_matrix = analytical_solution_KMB_matrix(t_size, x_size, dt, h);

        if (error_case) {
            auto kmb_error = calculate_absolute_error(KMB_matrix, KMB_analytical_matrix, t_size, x_size, t_size - 1);
            save(file, log10(h), log10(kmb_error));
        }
        if (analytical_case) {
            int t1 = t_size / 4;
            double time_1 = t1 * dt;
            int t2 = t_size / 2;
            double time_2 = t2 * dt;
            int t3 = t_size - 1;
            double time_3 = t3 * dt;
            file << "x \t ";
            file << "analitycal w t = " << time_1 << "\t";
            file << "KMB w t = " << time_1 << "\t";
            file << "analitycal w t = " << time_2 << "\t";
            file << "KMB w t = " << time_2 << "\t";
            file << "analitycal w t = " << time_3 << "\t";
            file << "KMB w t = " << time_3 << "\n";
            double x_iteration = r;
            for (int x = 1; x < x_size; x++) {
//                file << x_iteration << "\t" << KMB_analytical_matrix[t1][x] << "\t" << KMB_matrix[t1][x] << "\t"
//                     << KMB_analytical_matrix[t2][x] << "\t" << KMB_matrix[t2][x] << "\t"
//                     << KMB_analytical_matrix[t3][x] << "\t" << KMB_matrix[t3][x] << "\n";
                auto error = calculate_absolute_error(KMB_matrix, KMB_analytical_matrix, t_size, x_size, x);
                file << x_iteration << "\t" << log10(error) << "\n";
                x_iteration += h;
            }
        }
        if (time_case) {
            double t_iterator = dt;
            file << "t \t " << "error \n";
            for (int t = 1; t < t_size; t++) {
                auto error = calculate_absolute_error(KMB_matrix, KMB_analytical_matrix, t_size, x_size, t);
                file << t_iterator << "\t" << log10(error) << "\n";
                t_iterator += dt;
            }
        }
    }

    delete[] KMB_matrix;
    delete[] KMB_analytical_matrix;

    file.close();
}

void PMCN_THOMAS(double start_h, double stop_h, double step, const string &file_path, bool error_case,
                 bool analytical_case, bool time_case) {
    double dt;
    ofstream file;
    file.open(file_path);
    double **PMCN_matrix;
    double **PMCN_analytical_matrix;

    for (double h = start_h; h >= stop_h; h -= step) {
        // counting section
        dt = (h * h) * PMCN_lambda;
        t_size = floor(t_max / dt);
        x_size = floor(a / h);
        PMCN_matrix = init_matrix(t_size, x_size);
        init_matrix_with_conditions(PMCN_matrix, t_size, x_size, dt, h);
        PMCN_method(PMCN_matrix, h, dt, t_size, x_size);
        PMCN_analytical_matrix = analytical_solution_KMB_matrix(t_size, x_size, dt, h);

        if (error_case) {
            auto pmcn_error = calculate_absolute_error(PMCN_matrix, PMCN_analytical_matrix, t_size, x_size, t_size - 1);
            save(file, log10(h), log10(pmcn_error));
        }
        if (analytical_case) {
            int t1 = t_size / 4;
            double time_1 = t1 * dt;
            int t2 = t_size / 2;
            double time_2 = t2 * dt;
            int t3 = t_size - 1;
            double time_3 = t3 * dt;
            file << "x \t ";
            file << "analitycal w t = " << time_1 << "\t";
            file << "PCMN_THOMAS w t = " << time_1 << "\t";
            file << "analitycal w t = " << time_2 << "\t";
            file << "PCMN_THOMAS w t = " << time_2 << "\t";
            file << "analitycal w t = " << time_3 << "\t";
            file << "PCMN_THOMAS w t = " << time_3 << "\n";
            double x_iteration = r;
            for (int x = 0; x < x_size; x++) {
                file << x_iteration << "\t" << PMCN_analytical_matrix[t1][x] << "\t" << PMCN_matrix[t1][x] << "\t"
                     << PMCN_analytical_matrix[t2][x] << "\t" << PMCN_matrix[t2][x] << "\t"
                     << PMCN_analytical_matrix[t3][x] << "\t" << PMCN_matrix[t3][x] << "\n";
                x_iteration += h;
            }
        }
        if (time_case) {
            double t_iterator = dt;
            file << "t \t " << "error \n";
            for (int t = 1; t < t_size; t++) {
                auto error = calculate_absolute_error(PMCN_matrix, PMCN_analytical_matrix, t_size, x_size, t);
                file << t_iterator << "\t" << log10(error) << "\n";
                t_iterator += dt;
            }
        }
    }

    file.close();
}

void PMCN_LU(double start_h, double stop_h, double step, const string &file_path, bool error_case,
             bool analytical_case, bool time_case) {
    double dt;
    ofstream file;
    file.open(file_path);
    double **PMCN_LU_matrix;
    double **PMCN_LU_analytical_matrix;

    for (double h = start_h; h >= stop_h; h -= step) {
        // calculation sectino
        dt = (h * h) * PMCN_lambda;
        t_size = floor(t_max / dt);
        x_size = floor(a / h);
        PMCN_LU_matrix = init_matrix(t_size, x_size);
        init_matrix_with_conditions(PMCN_LU_matrix, t_size, x_size, dt, h);
        PMCN_LU_method(PMCN_LU_matrix, h, dt, t_size, x_size);
        PMCN_LU_analytical_matrix = analytical_solution_KMB_matrix(t_size, x_size, dt, h);

        if (error_case) {
            auto pmcn_lu_error = calculate_absolute_error(PMCN_LU_matrix, PMCN_LU_analytical_matrix, t_size, x_size,
                                                          t_size - 1);
            save(file, log10(h), log10(pmcn_lu_error));
        }
        if (analytical_case) {
            int t1 = t_size / 4;
            double time_1 = t1 * dt;
            int t2 = t_size / 2;
            double time_2 = t2 * dt;
            int t3 = t_size - 1;
            double time_3 = t3 * dt;
            file << "x \t ";
            file << "analitycal w t = " << time_1 << "\t";
            file << "PCMN_LU w t = " << time_1 << "\t";
            file << "analitycal w t = " << time_2 << "\t";
            file << "PCMN_LU w t = " << time_2 << "\t";
            file << "analitycal w t = " << time_3 << "\t";
            file << "PCMN_LU w t = " << time_3 << "\n";
            double x_iteration = r;
            for (int x = 0; x < x_size; x++) {
                file << x_iteration << "\t" << PMCN_LU_analytical_matrix[t1][x] << "\t" << PMCN_LU_matrix[t1][x] << "\t"
                     << PMCN_LU_analytical_matrix[t2][x] << "\t" << PMCN_LU_matrix[t2][x] << "\t"
                     << PMCN_LU_analytical_matrix[t3][x] << "\t" << PMCN_LU_matrix[t3][x] << "\n";
                x_iteration += h;
            }
        }
        if (time_case) {
            double t_iterator = dt;
            file << "t \t " << "error \n";
            for (int t = 1; t < t_size; t++) {
                auto error = calculate_absolute_error(PMCN_LU_matrix, PMCN_LU_analytical_matrix, t_size, x_size, t);
                file << t_iterator << "\t" << log10(error) << "\n";
                t_iterator += dt;
            }
        }


    }

    file.close();
}


int main() {
    //błedy KMB w odniesiu do h
    //KMB(0.5, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\KMB\KMB_error.txt)", true, false, false);

//     //rozwiązania analityczne i KMB dla dobrej dokładności
    KMB(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\KMB\tmp.txt)", false, true, false);
//
//     //błedy KMB w odniesieniu do t
//    KMB(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\KMB\KMB_error_time.txt)", false, false, true);
//
//
//
//
//     //błedy PMCN THOMAS w odniesieniu do h
//    PMCN_THOMAS(0.5, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN\PMCN_THOMAS_error.txt)", true, false, false);
//
//     //rozwiązania analityczne i PMCN THOMS dla dobrej dokładności
//    PMCN_THOMAS(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN\PMCN_THOMAS_functions.txt)", false, true, false);
//
//     //błedy PMCN THOMAS w odniesieniu do t
//    PMCN_THOMAS(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN\PMCN_THOMAS_error_time.txt)", false, false, true);
//
//
//
//
//     //błedy PMCN LU w odniesiu do h
//    PMCN_LU(0.5, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN_LU\PMCN_LU_error.txt)", true, false, false);
//
//     //rozwiązania analityczne i PMCN LU dla dobrej dokładności
//    PMCN_LU(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN_LU\PMCN_LU_functions.txt)", false, true, false);
//
//     //błedy PMCN LU w odniesieniu do t
//   PMCN_LU(0.01, 0.01, 0.01, R"(C:\studia\sem4\MO\MO_lab1_2\lab_11\PMCN_LU\PMCN_LU_error_time.txt)", false, false, true);

}

