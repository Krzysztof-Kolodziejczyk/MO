#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "Thomas.h"

using namespace std;

double alfa = 0., beta = 1., gamma = -1.;
double px = 1., qx = 0., rx = -4;
double fi = 0., psi = 1., theta = 0.;
double startX = 0., endX = 1.;
double functionDataThreshold = 110;

double sx(double x) {
    return -1. * x;
}

double U(double x) {
    double nominator = exp(2. - 2. * x) - 4. * exp(4. - 2. * x) + 4. * exp(2. * x) - exp(2. + 2. * x) - x + x * exp(4.);
    double denominator = 4. - 4. * exp(4.);
    return nominator / denominator;
}

double **initTriDiagonalMatrix(int n) {
    auto triDiagonal = new double *[n];
    for (int i = 0; i < n; i++) {
        triDiagonal[i] = new double[3];
    }
    return triDiagonal;
}

double maxFromVector(double *vector, int n) {
    double max = 0.;
    for (int i = 0; i < n; i++) {
        if (fabs(vector[i]) > max) {
            max = vector[i];
        }
    }
    return max;
}

void sendFunctionDataToFile(double xi, double approx, double analytic, ofstream &file) {
    file << setprecision(16) << xi << "\t" << approx << "\t" << analytic << endl;
}

// n-liczba podziałów przedziału [0,1]
double conventional(int n) {
    double h = (endX - startX) / (n - 1);

    auto errors = new double[n];

    auto analyticRes = new double[n];
    analyticRes[0] = U(startX);
    analyticRes[n - 1] = U(endX);

    auto bVector = new double[n];
    bVector[0] = -1. * gamma;
    bVector[n - 1] = -1. * theta;

    auto triDiagonal = initTriDiagonalMatrix(n);
    // warunek brzegowy w x = 0
    triDiagonal[0][0] = NULL;
    triDiagonal[0][1] = beta - alfa / h;
    triDiagonal[0][2] = alfa / h;
    // warunek brzegowy w x = 1
    triDiagonal[n - 1][0] = -1. * fi / h;
    triDiagonal[n - 1][1] = fi / h + psi;
    triDiagonal[n - 1][2] = NULL;
    // warunki wynikające z równania
    double xi = startX;
    for (int i = 1; i < n - 1; i++) {
        xi += h;
        triDiagonal[i][0] = px / (h * h) - qx / (2. * h);
        triDiagonal[i][1] = rx - 2. * px / (h * h);
        triDiagonal[i][2] = px / (h * h) + qx / (2. * h);
        bVector[i] = -1. * sx(xi);
        analyticRes[i] = U(xi);
    }

    auto thomasRes = Thomas::thomas(triDiagonal, n);
    Thomas::vectorB(thomasRes, bVector, n);
    auto xRes = Thomas::solve(thomasRes, bVector, n);

    if(n == functionDataThreshold){
        ofstream fileCon;
        fileCon.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_9\conFun.txt)");
        xi = startX;
        for(int i=0; i<n; i++){
            sendFunctionDataToFile(xi, xRes[i], analyticRes[i], fileCon);
            xi += h;
        }
    }

    for (int i = 0; i < n; i++) {
        errors[i] = xRes[i] - analyticRes[i];
    }

    delete[] triDiagonal;
    delete[] bVector;
    delete[] thomasRes;

    return log10(fabs(maxFromVector(errors, n)));

}

// n-liczba podziałów przedziału [0,1]
double numerow(int n) {

    auto errors = new double[n];
    double h = (endX - startX) / (n - 1);

    auto analyticRes = new double[n];
    analyticRes[0] = U(startX);
    analyticRes[n - 1] = U(endX);

    auto bVector = new double[n];
    bVector[0] = -1. * gamma;
    bVector[n - 1] = -1. * theta;

    auto triDiagonal = initTriDiagonalMatrix(n);
    // warunek brzegowy w x = 0
    triDiagonal[0][0] = NULL;
    triDiagonal[0][1] = beta - alfa / h;
    triDiagonal[0][2] = alfa / h;
    // warunek brzegowy w x = 1
    triDiagonal[n - 1][0] = -1. * fi / h;
    triDiagonal[n - 1][1] = fi / h + psi;
    triDiagonal[n - 1][2] = NULL;

    // warunki wynikające z równania
    double xi = startX;
    for (int i = 1; i < n - 1; i++) {
        xi += h;
        triDiagonal[i][0] = px / (h * h) + rx / 12.;
        triDiagonal[i][1] = (-2.0 * px) / (h * h) + rx * (10. / 12.);
        triDiagonal[i][2] = px / (h * h) + rx / 12.;
        bVector[i] = -1. * sx(xi - h) / 12. + -1. * (10. / 12.) * sx(xi) + -1. * sx(xi + h) / 12.;
        analyticRes[i] = U(xi);
    }

    auto thomasRes = Thomas::thomas(triDiagonal, n);
    Thomas::vectorB(thomasRes, bVector, n);
    auto xRes = Thomas::solve(thomasRes, bVector, n);

    if(n == functionDataThreshold){
        ofstream fileCon;
        fileCon.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_9\numFun.txt)");
        xi = startX;
        for(int i=0; i<n; i++){
            sendFunctionDataToFile(xi, xRes[i], analyticRes[i], fileCon);
            xi += h;
        }
    }

    for (int i = 0; i < n; i++) {
        errors[i] = xRes[i] - analyticRes[i];
    }

    delete[] triDiagonal;
    delete[] bVector;
    delete[] thomasRes;

    return log10(fabs(maxFromVector(errors, n)));
}

int main() {
    ofstream fileCon, fileNum;
    fileCon.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_9\conErr.txt)");
    fileNum.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_9\numErr.txt)");
    for (int i = 10; i < 30000; i += 50) {
        fileCon << log10((endX - startX) / (i - 1)) << "\t" << conventional(i) << endl;
        fileNum << log10((endX - startX) / (i - 1)) << "\t" << numerow(i) << endl;
    }
    return 0;
}