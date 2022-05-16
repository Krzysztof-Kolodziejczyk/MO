#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

#define M_PI		3.14159265358979323846

int ITERATIONS = 100;

template<typename T>
T function(T x) {
    return sin(x);
}

template<typename T>
T derivativeFun(T x) {
    return cos(x);
}

// may be used for forward, backward and central difference
template<typename T>
T twoPointsDiff(T x1, T x2, T step) {
    return (function(x2) - function(x1)) / step;
}

template<typename T>
T threePointsForwardDiff(T x0, T x1, T x2, T step) {
    return (T(-3) * function(x0) + T(4) * function(x1) - function(x2)) / (T(2) * step);
}

template<typename T>
T threePointsBackwardDiff(T x0, T x1, T x2, T step) {
    return (function(x0) + T(-4) * function(x1) + T(3) * function(x2)) / (T(2) * step);
}

template<typename T>
void saveResults(T *results, int n, ofstream &file){
    file << setprecision(16);
    for(int i=0; i<n; i++){
        file <<  setw(20) <<  results[i] << "\t";
    }
    file << "\n";
}

template<typename T>
void calculate() {
    T left = T(0);
    T right = T(M_PI) / 2.;
    T center = T(M_PI) / 4.;
    T h = T(0.1);
    auto results = new T[10];
    ofstream file;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_8\data.dat)");
    file <<  "#step\t0 forward (2)\t0 forward (3)\tpi/4 central (2)\tpi/4 forward (2)\tpi/4 backward (2)\tpi/4 forward (3)\tpi/4 backward (3)\tpi/2 backward (2)\tpi/2 backward (3)\n";
    for(int i=0; i<ITERATIONS; i++){
        results[0] = log10(h);
        results[1] = log10(fabs(derivativeFun(left) - twoPointsDiff(left, left + h, h))); // forward two points
        results[2] = log10(fabs(derivativeFun(left) - threePointsForwardDiff(left, left + h, left + T(2.0) * h, h))); // forward three points
        results[3] = log10(fabs(derivativeFun(center) - twoPointsDiff(center - h, center + h, T(2.0) * h))); // central two points
        results[4] = log10(fabs(derivativeFun(center) - twoPointsDiff(center, center + h, h))); // forward two points
        results[5] = log10(fabs(derivativeFun(center) - twoPointsDiff(center -h, center, h))); // backward two points
        results[6] = log10(fabs(derivativeFun(center) - threePointsForwardDiff(center, center + h, center + T(2.0) * h, h))); // forward three points
        results[7] = log10(fabs(derivativeFun(center) - threePointsBackwardDiff(center - T(2.0) * h, center - h, center, h))); // backward three points
        results[8] = log10(fabs(derivativeFun(right) - twoPointsDiff(right -h , right, h))); // backward two points
        results[9] = log10(fabs(derivativeFun(right) - threePointsBackwardDiff(right - T(2.0) * h, right -h, right, h))); // backward three points

        saveResults(results, 10, file);
        h *= 0.5;
    }
    file.close();
}

int main() {
    calculate<double>();
    return 0;
}