#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

/*
 *  Zaimplementuj w języku „C/C++” algorytm obliczający przybliżone wartości funkcji
    f(x) = [1 - exp(-x)]/x dla x E [10^30, 10^9], korzystając z funkcji standardowej exp(). W oparciu o
    zbiór dokładnych wartości tej funkcji, udostępniony przez prowadzącego zajęcia, zbadaj jak
    zmieniają się błędy względne przybliżenia funkcji w tym algorytmie, w zależności od x. W tym celu
    wykonaj rysunek przedstawiający zależność logarytmu dziesiętnego z bezwzględnej wartości błędu
    względnego od logarytmu dziesiętnego z x. Z wykresu odczytaj zakres zmiennej x, w którym błąd
    względny pozostaje na poziomie błędu reprezentacji, oraz zakres zmiennej x, w którym błąd
    względny jest większy. Wyjaśnij przyczynę obserwowanych zmian błędów. Na tej podstawie
    zaproponuj alternatywny sposób obliczania wartości funkcji f(x) w sytuacjach gdy obserwowany
    błąd jest duży. Dokonaj stosownej modyfikacji programu, tak aby uzyskać błąd względny na
    poziomie błędu reprezentacji (czyli tzw. dokładność maszynową) dla dowolnego x E [10^30, 10^9].
    W obliczeniach zastosuj zmienne podwójnej precyzji. Do wykonania rysunku w tym ćwiczeniu (a
    także w niektórych dalszych ćwiczeniach) najlepiej użyć programu GNUPLOT (dostępnego za
    darmo z Internetu).
 */

const int X_NBR = 3900;
double X_array[X_NBR];
double log_array[X_NBR];
double exp_array[X_NBR];


void build_input_arrays() {
    ifstream file;
    string line;
    file.open("C:\\studia\\MO_lab1_2\\zad_2\\oneminexp_ref.txt");

    int index = 0;
    double x;
    while (getline(file, line)) {
        istringstream iss(line);
        string nbr;
        int i = 0;
        while (iss >> nbr) {
            try {
                x = stod(nbr, nullptr);
                if (i == 0) {
                    log_array[index] = x;
                } else if (i == 1) {
                    X_array[index] = x;
                } else if (i == 2) {
                    exp_array[index] = x;
                    index++;
                }
                i++;
            } catch (exception ex) {
                break;
            }
        }
    }
}

double apprioximate_fun(double x) {
    return (1.0 - exp(-1.0 * x)) / x;
}

double relative_error(double af, double ef) {
    return (abs(af - ef) / ef);
}

double get_epsilon_from_double() {
    double e = 1.0;
    while (e / 2.0 + 1.0 > 1.0) {
        e = e / 2.0;
    }
    return e;
}

// w pliku log_outputs_expansion czasem wystepuje -inf. dzieje sie tak dlatego
// że roznica f_obliczone i f_prawidlowe wynosi zero (takie same bity) co skutkuje
// log10(0) = -inf. niekiedy wystepuje tez mala liczba rzedu 10^-16 ktora moze byc
// zwiazana z reprezentacja zmiennoprzecinkowa liczby (arytmetyka ni)
double fixed_function(double x) {
    double sum = 0.;
    double factor = 1.;
    double k = 2.0;
    double e = get_epsilon_from_double();

    while (abs(factor) > e) {
        sum += factor;
        factor *= (-1.0 * (x) / k);
        k += 1.0;
    }
    return sum;
}

int main() {

    build_input_arrays();

    ofstream approximate_result;
    approximate_result.open("C:\\studia\\MO_lab1_2\\zad_2\\approximate_result.txt", ofstream::out);
    ofstream fixed_result;
    fixed_result.open("C:\\studia\\MO_lab1_2\\zad_2\\fixed_result.txt", ofstream::out);
    ofstream relative_erros_approx;
    relative_erros_approx.open("C:\\studia\\MO_lab1_2\\zad_2\\relative_erros_approx.txt", ofstream::out);
    ofstream relative_erros_fixed;
    relative_erros_fixed.open("C:\\studia\\MO_lab1_2\\zad_2\\relative_erros_fixed.txt", ofstream::out);

    double ar, fr;
    for (int i = 0; i < X_NBR; i++) {

        // result
        ar = apprioximate_fun(X_array[i]);
        approximate_result << setprecision(20) << ar << endl;

        // fixed result
        if (X_array[i] <= 1.) {
            fr = fixed_function(X_array[i]);
        } else {
            fr = apprioximate_fun(X_array[i]);
        }
        //r = fixed_function(X_array[i]);
        fixed_result << setprecision(20) << fr << endl;

        // relative errors from approx
        relative_erros_approx << setprecision(20) << log_array[i] << " " << log10(abs(relative_error(ar, exp_array[i])))
                              << endl;

        // relative errors from fixed
        relative_erros_fixed << setprecision(20) << log_array[i] << " " << log10(abs(relative_error(fr, exp_array[i])))
                             << endl;
    }
}