#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

const int X_NBR = 3900;
double X_array[X_NBR];
double log_array[X_NBR];
double exp_array[X_NBR];


void build_input_arrays(){
    ifstream file;
    string line;
    file.open("C:\\studia\\MO_lab1_2\\zad_2\\oneminexp_ref.txt");

    int index = 0;
    double x;
    while(getline(file, line)){
        istringstream iss(line);
        string nbr;
        int i = 0;
        while(iss >> nbr) {
            try{
                x = stod(nbr, nullptr);
                if(i == 0){
                    log_array[index] = x;
                }
                else if (i == 1){
                    X_array[index] = x;
                }
                else if (i == 2){
                    exp_array[index] = x;
                    index++;
                }
                i++;
            }catch(exception ex){
                break;
            }
        }
    }
}

double apprioximate_fun(double x){
    return (1. - exp(-1. * x))/x;
}

double relative_error(double af, double ef){
    return (abs(af - ef)/ ef);
}

double get_epsilon_from_double(){
    double e = 1.;
    while(1.0 + e/2.0 > 1.0){
        e = e/2;
    }
    return e;
}

double fixed_function(double x){
    double result = 0.;
    double k = 1.;
    int i = 1;
    double e = get_epsilon_from_double();

    while( abs(k) > e){
        result += k;
        k *= (-1.0 * (x) / ((double)i + 1.0));
        i++;
    }
    return result;
}

int main(){

    build_input_arrays();

    ofstream approximate_result;
    approximate_result.open ("C:\\studia\\MO_lab1_2\\zad_2\\approximate_result.txt", ofstream::out);
    ofstream fixed_result;
    fixed_result.open ("C:\\studia\\MO_lab1_2\\zad_2\\fixed_result.txt", ofstream::out);
    ofstream relative_erros_approx;
    relative_erros_approx.open ("C:\\studia\\MO_lab1_2\\zad_2\\relative_erros_approx.txt", ofstream::out);
    ofstream relative_erros_fixed;
    relative_erros_fixed.open ("C:\\studia\\MO_lab1_2\\zad_2\\relative_erros_fixed.txt", ofstream::out);

    double ar, fr;
    for(int i=0; i<X_NBR; i++){

        // result
        ar = apprioximate_fun(X_array[i]);
        approximate_result << setprecision(20) << ar << endl;

        // fixed result
        if(X_array[i] <= 1.){
            fr = fixed_function(X_array[i]);
        }else{
            fr = apprioximate_fun(X_array[i]);
        }
        //r = fixed_function(X_array[i]);
        fixed_result << setprecision(20) << fr << endl;

        // relative errors from approx
        relative_erros_approx << setprecision(20) << log_array[i] << " " <<  log10(abs(relative_error(ar, exp_array[i]))) << endl;

        // relative errors from fixed
        relative_erros_fixed << setprecision(20) << log_array[i] << " " << log10(abs(relative_error(fr, exp_array[i]))) << endl;
    }
}