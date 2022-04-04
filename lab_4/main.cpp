#include <iostream>
#include <cmath>
#include<iomanip>
#include<algorithm>

using namespace std;

double NMAX = 200;
double TOLX = 10e-15;
double TOLF = 10e-15;

double **jacobiMatrix;
double **tmpMatrix;
double functionVector[3];
double xn[3];
double delta[3];

double det, det_x, det_y, det_z;

// funkcja 1
double f1(double x, double y, double z) {
    return x * x + y * y + z * z - 2;
}

double f1_dx(double x) {
    return 2.0 * x;
}

double f1_dy(double y) {
    return 2.0 * y;
}

double f1_dz(double z) {
    return 2.0 * z;
}

// funkcja 2
double f2(double x, double y) {
    return x * x + y * y - 1;
}

double f2_dx(double x) {
    return 2.0 * x;
}

double f2_dy(double y) {
    return 2.0 * y;
}

double f2_dz() {
    return 0;
}

// funkcja 3
double f3(double x, double y) {
    return x * x - y;
}

double f3_dx(double x) {
    return 2.0 * x;
}

double f3_dy() {
    return -1.0;
}

double f3_dz() {
    return 0;
}

double max(double wektor[3]) {
    return (max(max(wektor[0], wektor[1]), wektor[2]));
}

void buildJacobiMatrix(double x, double y, double z) {
    jacobiMatrix[0][0] = f1_dx(x);
    jacobiMatrix[0][1] = f1_dy(y);
    jacobiMatrix[0][2] = f1_dz(z);

    jacobiMatrix[1][0] = f2_dx(x);
    jacobiMatrix[1][1] = f2_dy(y);
    jacobiMatrix[1][2] = f2_dz();

    jacobiMatrix[2][0] = f3_dx(z);
    jacobiMatrix[2][1] = f3_dy();
    jacobiMatrix[2][2] = f3_dz();
}

void buildFunctionVector(double x, double y, double z) {
    functionVector[0] = f1(x, y, z);
    functionVector[1] = f2(x, y);
    functionVector[2] = f3(x, y);
}


double det3x3(double **matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[0][2]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

double **init3x3matrix() {
    auto **res = new double *[3];
    for (int i = 0; i < 3; i++) {
        res[i] = new double[3];
    }
    return res;
}

void resetTmpMatrix() {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            tmpMatrix[i][j] = jacobiMatrix[i][j];
        }
    }
}

void fillTmpMatrixColumnWithFunctionVector(int column) {
    resetTmpMatrix();
    for (int i = 0; i < 3; i++) {
        tmpMatrix[i][column] = functionVector[i];
    }
}

void calculateDelta(){
    delta[0] = det_x / det;
    delta[1] = det_y / det;
    delta[2] = det_z / det;
}

void printVectorXn() {
    cout << "[" << setprecision(7) << xn[0] << ", " << xn[1] << " ," << xn[2] << "]" << endl;
}

int main() {
    xn[0] = 1.0, xn[1] = 1.0, xn[2] = 1.0;
    int i = 0;

    jacobiMatrix = init3x3matrix();
    tmpMatrix = init3x3matrix();
    do {
        buildJacobiMatrix(xn[0], xn[1], xn[2]);
        buildFunctionVector(xn[0], xn[1], xn[2]);
        det = det3x3(jacobiMatrix);
        fillTmpMatrixColumnWithFunctionVector(0);
        det_x = det3x3(tmpMatrix);
        fillTmpMatrixColumnWithFunctionVector(1);
        det_y = det3x3(tmpMatrix);
        fillTmpMatrixColumnWithFunctionVector(2);
        det_z = det3x3(tmpMatrix);
        calculateDelta();

        xn[0] -= delta[0];
        xn[1] -= delta[1];
        xn[2] -= delta[2];

        cout << "i = " << i << "\t" << setprecision(7) << "x_error = " << max(delta) << "\t" << "reziduum = " << max(functionVector) << "\t";
        printVectorXn();
        i++;
    } while (((max(delta) > TOLX) || (max(functionVector) > TOLF)) && i < NMAX);
}