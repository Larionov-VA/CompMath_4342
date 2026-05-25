#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#include <vector>
#include <string>

using namespace std;

double f(double x) {
    return cos(x * x * x);
}

double rectangle(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + (i + 0.5) * h;
        sum += f(x);
    }
    return h * sum;
}

double trapezoid(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }
    return h * sum;
}

double simpson(double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += (i % 2 == 1) ? 4.0 * f(x) : 2.0 * f(x);
    }
    return (h / 3.0) * sum;
}

struct Result {
    string method;
    double eps;
    int n;
    double h;
    double integral;
    double rungeIntegral;
    double error;
    int iterations;
};

Result rungeRule(const string& methodName,
                 function<double(double,double,int)> formula,
                 int order, double a, double b, double eps) {
    int n = 2;
    double oldIntegral = formula(a, b, n);
    int iter = 0;
    while (true) {
        n *= 2;
        iter++;
        double newIntegral = formula(a, b, n);
        double denom = pow(2.0, order) - 1.0;
        double error = fabs(newIntegral - oldIntegral) / denom;
        double rungeVal = newIntegral + (newIntegral - oldIntegral) / denom;
        if (error <= eps) {
            return {methodName, eps, n, (b-a)/n, newIntegral, rungeVal, error, iter};
        }
        oldIntegral = newIntegral;
    }
}

void printResult(const Result& r) {
    cout << fixed << setprecision(10);
    cout << "Метод: " << r.method << "\n\n";
    cout << "n = " << r.n << "\n";
    cout << "h = " << r.h << "\n";
    cout << "I по формуле = " << r.integral << "\n";
    cout << "I с уточнением по Рунге = " << r.rungeIntegral << "\n";
    cout << "Оценка погрешности = " << r.error << "\n";
    cout << "Удвоений n = " << r.iterations << "\n";
    cout << "--------------------------------------------\n";
}

int main() {
    double a = 0.0, b = 1.0;
    vector<double> epsilons = {0.01, 0.001, 0.0001};

    for (double eps : epsilons) {
        cout << "\nТочность eps = " << eps << "\n\n";
        printResult(rungeRule("Прямоугольники", rectangle, 2, a, b, eps));
        printResult(rungeRule("Трапеции", trapezoid, 2, a, b, eps));
        printResult(rungeRule("Симпсон", simpson, 4, a, b, eps));
    }
    return 0;
}