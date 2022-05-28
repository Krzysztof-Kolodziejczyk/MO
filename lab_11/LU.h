

#ifndef MO_LAB1_2_LU_H
#define MO_LAB1_2_LU_H


class LU {

public:
    static int *gauss(double **A, int n);
    static void vectorB(double *b, double **A, const int *indexes, int n);
    static double *solve(double **A, const double *b, const int *indexes, int n);
    static void gauss_recur(double **A, int *indexes, int n, int p = 0);
};


#endif //MO_LAB1_2_LU_H
