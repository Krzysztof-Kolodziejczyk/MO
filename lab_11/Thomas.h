

#ifndef MO_LAB1_2_THOMAS_H
#define MO_LAB1_2_THOMAS_H


class Thomas {

public:

    static double **thomas(double **A, int n);

    static void vectorB(double **matrixThomas, double *b, int n);

    static double *solve(double **matrixThomas, const double *b, int n);

};


#endif //MO_LAB1_2_THOMAS_H
