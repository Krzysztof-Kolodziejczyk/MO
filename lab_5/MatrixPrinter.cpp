

#include "MatrixPrinter.h"

using namespace std;

void MatrixPrinter::printMatrix(double **m, int n) {
    for(int j=0; j<49; j++){
        //cout << "=";
    }
    //cout << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout  << setw(8) << ceil(m[i][j] * 1000.0) / 1000.0 << "   ";
        }
//        cout << "|" << endl;
//        for(int j=0; j<49; j++){
//            cout << "=";
//        }
        cout << endl;
    }
    cout << endl << endl;
}
