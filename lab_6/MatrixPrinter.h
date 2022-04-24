

#ifndef MO_LAB1_2_MATRIXPRINTER_H
#define MO_LAB1_2_MATRIXPRINTER_H

#include "iostream"
#include "iomanip"
#include "string"

using namespace std;

class MatrixPrinter {

    static void printSeparator(int n);
public:
    static void printMatrix(double **A, const int *indexes, int n, int m, const string& name);

    static void printVector(const double *v, const int *indexes, int n, const string& name);
};


#endif //MO_LAB1_2_MATRIXPRINTER_H
