#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <fstream>

using namespace std;

double STOP = 1.;
double ITERS = 1000.;

double fAnalytic(double t) {
    return 1. - exp(-10. * (t + atan(t)));
}

double BME() {
    double dt = STOP / ITERS;
    double analytic_result, t = 0, maxError = 0, error;
    double y = 0.;
    ofstream file;
    //file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_10\BME_stable.txt)");

    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_10\BME_unstable.txt)");

    while(t < STOP){
        analytic_result = fAnalytic(t+dt);
        y = y -1.* (((10. * t * t) + 20.) / (t * t + 1)) * (y - 1) * dt;
        t += dt;
        file << t << "\t" << y << "\t" << analytic_result << endl;
        error = fabs(y - analytic_result);
        if (error > maxError)
            maxError = error;
    }
    return error;
}

double PME() {
    double analytic_result, tk1, maxError = 0, error, nominator, denominator;
    double y = 0.;
    double dt = STOP / ITERS;
    tk1 = dt;

    //ofstream file;
    //file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_10\PME.txt)");

    while(tk1 < STOP){
        analytic_result = fAnalytic(tk1);
        nominator = dt * (10. * tk1 * tk1 + 20.) + y * (tk1 * tk1 + 1.);
        denominator = dt * (10. * tk1 * tk1 + 20.) + (tk1 * tk1 + 1.);
        y = nominator / denominator;
        //file << tk1 << "\t" << y << "\t" << analytic_result << endl;
        tk1 += dt;
        error = fabs(y - analytic_result);
        if (error > maxError)
            maxError = error;
    }
    return error;
}

double PMT() {
    double analytic_result, tk1, t=0, maxError = 0, error, nominator, denominator;
    double y = 0.;
    double dt = STOP / ITERS;
    tk1 = dt;
    double p0, p1;

    //ofstream file;
    //file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_10\PMT.txt)");

    while(t < STOP){
        analytic_result = fAnalytic(tk1);
        p0 = (10. * t * t + 20.) / (t * t + 1.);
        p1 = (10. * tk1 * tk1 + 20.) / (tk1 * tk1 + 1.);
        nominator = 2. * y / dt - p0 * (y - 1.) + p1;
        denominator = 2. / dt + p1;
        y = nominator / denominator;
        //file << tk1 << "\t" << y << "\t" << analytic_result << endl;
        t += dt;
        tk1 += dt;
        error = fabs(y - analytic_result);
        if (error > maxError)
            maxError = error;
    }
    return error;
}


int main() {
    ITERS = 10;
    ofstream file;
    file.open(R"(C:\studia\sem4\MO\MO_lab1_2\lab_10\Errors.txt)");
    for(ITERS=10; ITERS<1.e6; ITERS+=1000){
        // errors
        file << log10(STOP / ITERS) << "\t" << log10(BME()) << "\t" << log10(PME()) << "\t" << log10(PMT()) << endl;

    }
    file.close();

//    cout << BME() << endl;
//    cout << PME() << endl;
//    cout << PMT() << endl;
}

