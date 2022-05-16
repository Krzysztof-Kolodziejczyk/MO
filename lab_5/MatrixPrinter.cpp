

#include "MatrixPrinter.h"

void MatrixPrinter::printMatrix(double **m, const int *indexes, int n, const string& name) {
    cout << name << endl;
    printSeparator(49);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(indexes == nullptr){
                cout << "|" << setw(8) << setprecision(5) << m[i][j] << "   ";
            }else{
                cout << "|" << setw(8) << setprecision(5) << m[indexes[i]][j] << "   ";
            }
        }
        cout << "|" << endl;
        printSeparator(49);
    }
    cout << endl;
}

void MatrixPrinter::printL(double **m, const int *indexes, int n, const string &name) {
    cout << name << endl;
    printSeparator(49);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i == j){
                cout << "|" << setw(8) << 1. << "   ";
            }else if(i < j){
                cout << "|" << setw(8) << 0. << "   ";
            }
            else{
                cout << "|" << setw(8) << setprecision(5) << m[indexes[i]][j] << "   ";
            }
        }
        cout << "|" << endl;
        printSeparator(49);
    }
    cout << endl;
}

void MatrixPrinter::printU(double **m, const int *indexes, int n, const string &name) {
    cout << name << endl;
    printSeparator(49);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i > j){
                cout << "|" << setw(8) << 0. << "   ";
            }
            else{
                cout << "|" << setw(8) << setprecision(5) << m[indexes[i]][j] << "   ";
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
        if(indexes != nullptr){
            cout << "|" << setw(8) << setprecision(5) << v[indexes[i]] << "   ";
        }else{
            cout << "|" << setw(8) << setprecision(5) << v[i] << "   ";
        }

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

