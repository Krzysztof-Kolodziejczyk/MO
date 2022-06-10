#include<iostream>
#include <iomanip>
#include <cfloat>
#include "cmath"
#include <fstream>

using namespace std;

double a = -1.;
double b = 1.;
int N = 20;
double PI = 3.14159265358979323;

double f(double x) {
    return 1. / (1. + 10. * x * x * x * x * x * x);
}

double *c_vector(const double *xi, const double *fxi, int n) {
    auto *tmp = new double[n];
    for (int i = 0; i < n; i++) {
        tmp[i] = fxi[i];
    }
    auto *c = new double[n];
    c[0] = tmp[0];
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp[j] = (tmp[j + 1] - tmp[j]) / (xi[j + i] - xi[j]);
        }
        c[i] = tmp[0];
    }
//    for (int i = 0; i < N; i++) {
//        cout << c[i] << endl;
//    }
    return c;
}

double netwon(const double *c, double x, const double *xi, int n) {
    double res = c[n - 1] * (x - xi[n - 2]) + c[n - 2];
    for (int i = n - 3; i >= 0; i--) {
        res *= (x - xi[i]);
        res += c[i];
    }
    return res;
}


int main() {
    ofstream file, file1, file2;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_12\function.txt)");
    file1.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_12\normal.txt)");
    file2.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_12\czebyszew.txt)");
    double dx = fabs(b - a) / N;
    double x = a;
    auto *xi = new double[N];
    auto *fxi = new double[N];

    for (int i = 0; i < N; i++) {
        xi[i] = x;
        fxi[i] = f(x);
        x += dx;
    }

    auto c = c_vector(xi, fxi, N);
    dx = fabs(b - a) / 100;
    x = a;
    for (int i = 0; i < 100; i++) {
        file1 << setprecision(20) << x << "\t" << netwon(c, x, xi, N) << endl;
        file << setprecision(20) << x << "\t" << f(x) << endl;
        x += dx;
    }

    // wezły czebyszewa
    for (int i = 0; i < N; i++) {
        xi[i] = (b + a) / 2. + ((b - a) / 2.) * cos(PI * (2. * i + 1) / (2. * N - 1 + 2));
        fxi[i] = f(xi[i]);
    }

    c = c_vector(xi, fxi, N);
    dx = fabs(b - a) / 100;
    x = a;
    for (int i = 0; i < 100; i++) {
        file2 << setprecision(20) << x << "\t" << netwon(c, x, xi, N) << endl;
        x += dx;
    }
}

