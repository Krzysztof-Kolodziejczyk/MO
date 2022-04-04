#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

int ITER = 50;
double TOLX = 10e-10, TOLF = 10e-10;

double fun_1(double x) {
    return sin(x / 4.0) * sin(x / 4.0) - x;
}

double fi_fun_1(double x) {
    return sin(x / 4.0) * sin(x / 4.0);
}

double fun_1_dx(double x) {
    return 0.5 * sin(x / 4.0) * cos(x / 4.0) - 1.0;
}

double fun_2(double x) {
    return tan(2.0 * x) - x - 1;
}

double fi_fun_2(double x) {
    return tan(2.0 * x) - 1;
}

double fun_2_dx(double x) {
    return (2.0 * 1 / (cos(2.0 * x) * cos(2.0 * x))) - 1;
}


void picard() {
    cout << "\n-----------------------  Picard  -----------------------\n\n";
    cout << "fun 1 " << endl;
    int i = 0;
    double x = 0.5, x0;

    do{
        x0 = x;
        x = fi_fun_1(x);
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs(x - x0) << "\t reziduum = " << setw(10) <<  fabs(fun_1(x)) << endl;
    }while (i++ <= ITER && !(fabs(x - x0) < TOLX && fabs(fun_1(x)) < TOLF));

    cout << endl << endl << "fun 2 " << endl;
    i = 0, x = 0.5;
    do{
        x0 = x;
        x = fi_fun_2(x);
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs(x - x0) << "\t reziduum = " << setw(10) <<  fabs(fun_2(x)) << endl;
    }while (i++ <= ITER && !(fabs(x - x0) < TOLX && fabs(fun_1(x)) < TOLF));

    cout << "\n-----------------------  Picard - koniec -----------------------\n\n\n\n";
}


void bisekcja() {
    cout << "\n-----------------------  Bisekcja  -----------------------\n\n";
    cout << "fun 1 " << endl;
    double a = -10.0, b = 20.0, x;
    int i = 0;
    do {
        x = (a + b) / 2;
        if (fun_1(a) < 0 && fun_1(x) > 0 || fun_1(a) > 0 && fun_1(x) < 0) {
            b = x;
        } else if (fun_1(b) > 0 && fun_1(x) < 0 || fun_1(b) < 0 && fun_1(x) > 0) {
            a = x;
        } else {
            cout << "brak rozwiazania w danym przedziale" << endl;
            break;
        }
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs((b - a) / 2.) << "\t reziduum = " << setw(10) <<  fabs(fun_1(x)) << endl;
    } while (i++ <= ITER && !(fabs((b - a) / 2.) <= TOLX && fabs(fun_1(x)) <= TOLF));


    cout << endl << endl << "fun 2 " << endl;
    a = -0.5, b = 0.7;
    i = 0;

    do {
        x = (a + b) / 2;
        if (fun_2(a) < 0 && fun_2(x) > 0 || fun_2(a) > 0 && fun_2(x) < 0) {
            b = x;
        } else if (fun_2(b) > 0 && fun_2(x) < 0 || fun_2(b) < 0 && fun_2(x) > 0) {
            a = x;
        } else {
            cout << "brak rozwiazania w danym przedziale" << endl;
            break;
        }
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs((b - a) / 2.) << "\t reziduum = " << setw(10) <<  fabs(fun_2(x)) << endl;
    } while (i++ <= ITER && !(fabs((b - a) / 2.) < TOLX && fabs(fun_2(x)) < TOLF));

    cout << "\n-----------------------  Bisekcja - koniec  -----------------------\n\n\n\n";
}


void newton() {
    cout << "\n-----------------------  Newton  -----------------------\n\n";
    int i = 0;
    double x = 0.5, x0;

    cout << "fun 1 " << endl;
    do{
        x0 = x;
        x = x - (fun_1(x) / fun_1_dx(x));
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs(x - x0) << "\t reziduum = " << setw(10) <<  fabs(fun_1(x)) << endl;
    }while (i++ <= ITER && !(fabs(x - x0) < TOLX && fabs(fun_1(x)) < TOLF));

    cout << endl << endl << "fun 2 " << endl;
    i=0, x=0.5;
    do{
        x0 = x;
        x = x - (fun_2(x) / fun_2_dx(x));
        cout << "i = " << i << "\t x = " << setw(10) << x << "\t x_error = " << setw(10) << fabs(x - x0) << "\t reziduum = " << setw(10) <<  fabs(fun_2(x)) << endl;
    }while (i++ <= ITER && !(fabs(x - x0) < TOLX && fabs(fun_2(x)) < TOLF));

    cout << "\n-----------------------  Newton - koniec  -----------------------\n\n\n\n";
}


void sieczne() {
    cout << "\n-----------------------  Sieczne  -----------------------\n\n";
    int i = 0;
    double x0, x1 = 3.0, x2 = 7.0;
    cout << "fun 1 " << endl;
    do{
        x0 = x1;
        x1 = x2;
        x2 = x1 - (fun_1(x1) * (x1 - x0)) / (fun_1(x1) - fun_1(x0));
        cout << "i = " << i << "\t x = " << setw(10) << x2 << "\t x_error = " << setw(10) << fabs(x2 - x1) << "\t reziduum = " << setw(10) <<  fabs(fun_1(x2)) << endl;
    }while (i++ <= ITER && !(fabs(x2 - x1) < TOLX && fabs(fun_1(x2)) < TOLF));

    x1 = .5, x2 = .7;
    cout << endl << endl << "fun 2 " << endl;
    do{
        x0 = x1;
        x1 = x2;
        x2 = x1 - (fun_2(x1) * (x1 - x0)) / (fun_2(x1) - fun_2(x0));
        cout << "i = " << i << "\t x = " << setw(10) << x2 << "\t x_error = " << setw(10) << fabs(x2 - x1) << "\t reziduum = " << setw(10) <<  fabs(fun_2(x2)) << endl;
    }while (i++ <= ITER && !(fabs(x2 - x1) < TOLX && fabs(fun_2(x2)) < TOLF));

    cout << "\n-----------------------  Sieczne - koniec  -----------------------\n\n\n\n";
}


int main() {
    picard();
    bisekcja();
    newton();
    sieczne();
}