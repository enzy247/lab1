#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

/*int main() {
    complex<double> z1(-1, 1);
    complex<double> z2(-2, -2);
    complex<double> z3 = pow(complex<double>(-1, 1), 1.0/3.0);
    complex<double> z4 = pow(complex<double>(-1, 1), 4);

    complex<double> sum = z1 + z2;
    cout << "Сумма: " << sum << endl;

    complex<double> diff = z1 - z2;
    cout << "Разность: " << diff << endl;

    complex<double> prod = z1 * z2;
    cout << "Произведение: " << prod << endl;

    complex<double> quot = z1 / z2;
    cout << "Частное: " << quot << endl;
    
    cout << "Четвертая степень: " << z4 << endl;

    cout << "Корень третьей степени: " << z3 << endl;
    return 0;
}
*/

double f(double x) {
    return log(x) + pow(x + 1, 3);
}

double df(double x) {
    return 1 / x + 3 * pow(x + 1, 2);
}

double bisect(double a, double b, double eps) {
    double c;
    while ((b - a) >= eps) {
        c = (a + b) / 2;
        if (f(c) == 0.0)
            break;
        else if (f(c) * f(a) < 0)
            b = c;
        else
            a = c;
    }
    return c;
}

double iterate(double x0, double eps) {
    double x1;
    while (true) {
        x1 = f(x0);
        x0 = x1;
        if (fabs(x1 - x0) < eps) 
            break;
    } 
    return x1;
    
}

double chord(double a, double b, double eps) {
    double x0 = a, x1;
    while (fabs(f(x1)) >= eps) {
        x1 = x0 - f(x0) * (b - x0) / (f(b) - f(x0));
        if (f(x1) == 0.0) break;
        else if (f(x1) * f(a) < 0) b = x1;
        else a = x1;
        x0 = x1;
    } 
    return x1;
}

double newton(double x0, double eps) {
    double x1;
    while (fabs(f(x1)) >= eps) {
        x1 = x0 - f(x0) / df(x0);
        x0 = x1;
    } 
    return x1;
}

int main() {
    double eps = 1e-6;
    double guess = 1.0;
    double a = 0.1, b = 2.0;

    cout << "Метод дихотомии: "  << bisect(a, b, eps) << endl;
    cout << "Метод итераций: " << iterate(guess, eps) << endl;
    cout << "Метод хорд: " << chord(a, b, eps) << endl;
    cout << "Метод Ньютона: " << newton(guess, eps) << endl;

    return 0;
}