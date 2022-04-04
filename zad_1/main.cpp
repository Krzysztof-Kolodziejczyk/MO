#include<iostream>
#include <iomanip>
#include<cfloat>
#include "math.h"

using namespace std;

/*
 *  Napisz program w języku „C/C++”, umożliwiający „doświadczalne” wyznaczenie liczby bitów
    mantysy oraz tzw. epsylona maszynowego, dla zmiennych typu float i double, tj. najmniejszej liczby
     takiej, że fl( + 1) > 1. Aby znaleźć odpowiedź na pytanie jak napisać taki program, zacznij od
    wyjaśnienia kwestii jaki jest związek  z precyzją arytmetyki.
 */


/*
 * Epsilon maszynowy – wartość określająca precyzję obliczeń numerycznych wykonywanych na liczbach zmiennoprzecinkowych.
 * Jest to największa liczba nieujemna, której dodanie do jedności daje wynik równy 1
 *
 * precyzja arytmetyki - maksymalny bląd względy reprezentacji liczby znormalizowanej
 */
void float_precission(){
    int t = 0;
    // e na początku = 1 bo to jest 2^0
    float e = 1.f;

    while(e/2.f + 1.f > 1.f){
        t++;
        // dzielimy bo finalnie checemu uzyskać wynik 2^(-t)
        e = e/2.f;
    }
    cout << "mantissa bit nbr: " << setprecision(20) << t << endl;
    cout << "epsilon : " << setprecision(20) << e << endl;
    cout << "epsilon from lib : " << setprecision(20) << FLT_EPSILON   << endl;
}

void double_precission(){
    int t = 0;
    double e = 1.0;

    while(e/2.0 + 1.0 > 1.0){
        t++;
        e = e/2.0;
    }
    cout << "mantissa bit nbr: " << setprecision(20) << t << endl;
    cout << "epsilon : " << setprecision(20) << e << endl;
    cout << "epsilon from lib : " << setprecision(20) << DBL_EPSILON  << endl;
}

int main(){
    cout << "float\n";
    float_precission();
    cout  << "\ndouble\n";
    double_precission();

    cout <<endl << setprecision(20) <<  1.f + (float) 1.1920928955078125e-07;
    cout <<endl << setprecision(20) <<  1.f + (float) 1.1920928955078126e-07/2.f << endl;


    cout << setprecision(20) << 1.0 - exp(-1.0 * 1) << endl;
    cout << setprecision(20) << 1.0 -  exp(-1.0 * 0.01) << endl;
    cout << setprecision(20) << 1.0 - exp(-1.0 * 0.0001) << endl;
    cout << setprecision(20) << 1.0 -  exp(-1.0 * 0.000001) << endl;
    cout << setprecision(20) << 1.0 - exp(-1.0 * 0.000000001) << endl;
    cout << setprecision(20) << 1.0 -  exp(-1.0 * 0.00000000000000001) << endl;
}

