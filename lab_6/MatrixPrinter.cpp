

#include "MatrixPrinter.h"

void MatrixPrinter::printMatrix(double **A, const int *indexes, int n, int m, const string& name) {
    cout << name << endl;
    //printSeparator(61);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if(indexes == nullptr){
                cout << "|" << setw(15) << setprecision(11) <<  A[i][j];
            }else{
                cout << "|" << setw(15) << setprecision(11) << A[indexes[i]][j];
            }
        }
        cout << "|" << endl;
        //printSeparator(61);
    }
    cout << endl;
}

void MatrixPrinter::printVector(const double *v, const int *indexes, int n, const string& name) {
    cout << name << endl;
    //printSeparator(73);
    for(int i=0; i<n; i++){
        if(indexes == nullptr){
            cout << "|" << setw(15) << setprecision(11) << v[i];
        } else{
            cout << "|" << setw(15) << setprecision(11) << v[indexes[i]];
        }
    }
    cout << "|" << endl;
    //printSeparator(73);
    cout << endl;
}

//void MatrixPrinter::printSeparator(int N) {
//    for (int j = 0; j < N; j++) {
//        cout << "=";
//    }
//    cout << endl;
//}
