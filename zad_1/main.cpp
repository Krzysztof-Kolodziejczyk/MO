#include<iostream>
#include <iomanip>
#include<cfloat>

using namespace std;

void float_precission(){
    int t = 0;
    float e = 1.f;

    while(e/2 + 1.f > 1.f){
        t++;
        e = e/2;
    }
    cout << "mantissa bit nbr: " << setprecision(20) << t << endl;
    cout << "epsilon : " << setprecision(20) << e << endl;
    cout << "epsilon from lib : " << setprecision(20) << FLT_EPSILON   << endl;
}

void double_precission(){
    int t = 0;
    double e = 1.0;

    while(e/2 + 1.0 > 1.0){
        t++;
        e = e/2;
    }
    cout << "mantissa bit nbr: " << setprecision(20) << t << endl;
    cout << "epsilon : " << setprecision(20) << e << endl;
    cout << "epsilon from lib : " << setprecision(20) << DBL_EPSILON  << endl;
}

int main(){
    float_precission();
    cout << endl << endl;
    double_precission();
}

