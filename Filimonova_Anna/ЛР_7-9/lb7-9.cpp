#include <cmath>
#include <cstdio>

double f(double x) {
    return std::cos(x + x * x * x);
}

double rectangle(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x_mid = a + (i + 0.5) * h;
        sum += func(x_mid);
    }
    return sum * h;
}

double trapezoid(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += (func(a + i * h) + func(a + (i + 1) * h));
    }
    return sum * 0.5 * h;
}

double simpson(double (*func)(double), double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; i += 2) {
        sum += (func(a + i * h) + 4.0 * func(a + (i + 1) * h) + func(a + (i + 2) * h));
    }
    return sum * h / 3.0;
}

void integrate_with_runge(
    double (*method)(double (*)(double), double, double, int),
    const char* method_name,
    double divisor,           
    double a, double b,
    double eps,
    int precision)
{
    int n = 2;
    double I_n = method(f, a, b, n);

    while (true) {
        int n2 = n * 2;
        double I_2n = method(f, a, b, n2);
        double error = std::fabs(I_2n - I_n) / divisor;   

        if (error < eps) {
            double I_final = I_2n + (I_2n - I_n) / divisor;   
            printf("%s: I = %.*f, n = %d, h = %.*f, погрешность = %.*f\n",
                   method_name, precision, I_final,
                   n2, precision, (b - a) / n2,
                   precision, error);
            break;
        }
        I_n = I_2n;
        n = n2;
    }
}

int main() {
    const double PI = std::acos(-1.0);
    double a = PI / 2.0;
    double b = PI;

    printf("Интеграл от pi/2 до pi функции cos(x + x^3) dx\n\n");

    for (int k = 2; k <= 10; ++k) {
        double eps = pow(10.0, -k);   
        int prec = k + 1;             

        printf("Точность eps = %g\n", eps);
        integrate_with_runge(rectangle, "Прямоугольники", 3.0,  a, b, eps, prec);
        integrate_with_runge(trapezoid,  "Трапеции",       3.0,  a, b, eps, prec);
        integrate_with_runge(simpson,    "Симпсон",        15.0, a, b, eps, prec);
        printf("\n");
    }
    return 0;
}