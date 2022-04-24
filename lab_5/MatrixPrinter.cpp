

#include "MatrixPrinter.h"

void MatrixPrinter::printMatrix(double **m, const int *indexes, int n, const string& name) {
    cout << name << endl;
    printSeparator(49);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(indexes == nullptr){
                cout << "|" << setw(8) << ceil(m[i][j] * 1000.0) / 1000.0 << "   ";
            }else{
                cout << "|" << setw(8) << ceil(m[indexes[i]][j] * 1000.0) / 1000.0 << "   ";
            }
        }
        cout << "|" << endl;
        printSeparator(49);
    }
    cout << endl;
}

void MatrixPrinter::printVector(const double *v, const int *indexes, int n, const string& name) {
    cout << name << endl;
    printSeparator(49);
    for(int i=0; i<n; i++){
        cout << "|" << setw(8) << ceil(v[indexes[i]] * 1000.0) / 1000.0 << "   ";
    }
    cout << "|" << endl;
    printSeparator(49);
    cout << endl;
}

void MatrixPrinter::printSeparator(int n) {
    for (int j = 0; j < n; j++) {
        cout << "=";
    }
    cout << endl;
}
