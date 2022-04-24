

#include "MatrixPrinter.h"

void MatrixPrinter::printMatrix(double **A, const int *indexes, int n, int m, const string& name) {
    cout << name << endl;
    printSeparator(61);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if(indexes == nullptr){
                cout << "|" << setw(16) << ceil(A[i][j] * 1000000000.0) / 1000000000.0 << "   ";
            }else{
                cout << "|" << setw(16) << ceil(A[indexes[i]][j] * 10000000000.0) / 1000000000.0 << "   ";
            }
        }
        cout << "|" << endl;
        printSeparator(61);
    }
    cout << endl;
}

void MatrixPrinter::printVector(const double *v, const int *indexes, int n, const string& name) {
    cout << name << endl;
    printSeparator(73);
    for(int i=0; i<n; i++){
        if(indexes == nullptr){
            cout << "|" << setw(8) << ceil(v[i] * 1000.0) / 1000.0 << "   ";
        } else{
            cout << "|" << setw(8) << ceil(v[indexes[i]] * 1000.0) / 1000.0 << "   ";
        }
    }
    cout << "|" << endl;
    printSeparator(73);
    cout << endl;
}

void MatrixPrinter::printSeparator(int n) {
    for (int j = 0; j < n; j++) {
        cout << "=";
    }
    cout << endl;
}
