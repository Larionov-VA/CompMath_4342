#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>

using namespace std;

double F(double x)
{
    return sin(pow(x, 4) + 2.0 * pow(x, 3) + x * x);
}

double rectangleFormula(double a, double b, int n)
{
    double h = (b - a) / n; // Вычисляем шаг
    double sum = 0.0;       // Переменная для накопления суммы значений функции

    for (int i = 0; i < n; i++)
    {
        // Выбираем x в середине каждого отрезка
        double x = a + (i + 0.5) * h;
        sum += F(x);
    }

    return h * sum; // Площадь (ширина * сумма высот)
}

double trapezoidFormula(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.5 * (F(a) + F(b));

    for (int i = 1; i < n; i++)
    {
        double x = a + i * h;   // x внутри отрезка [a, b]
        sum += F(x);
    }

    return h * sum;
}

double simpsonFormula(double a, double b, int n)
{
    double h = (b - a) / n;
    // (h/3) * (f(a) + f(b) + 4*сумма_нч + 2*сумма_ч)
    double sum = F(a) + F(b);   // Значения на краях отрезка

    for (int i = 1; i < n; i++)
    {
        double x = a + i * h;

        if (i % 2 == 0)
            sum += 2.0 * F(x);  // Чет узлы имеют коэффициент 2
        else
            sum += 4.0 * F(x);
    }

    return h * sum / 3.0;
}

struct Result
{
    int n;
    double h;
    double valueH;
    double valueHalfH;
    double rungeCorrection;
    double integralValue;
};

typedef double (*Method)(double, double, int);

Result solveByRunge(Method method, double divider, int sign, double a, double b, double eps)
{
    int n = 2;

    while (true)
    {
        double valueH = method(a, b, n);
        double valueHalfH = method(a, b, 2 * n);

        double rungeCorrection = fabs(valueHalfH - valueH) / divider;

        double integralValue = valueHalfH + sign * rungeCorrection;

        if (rungeCorrection <= eps)
        {
            Result result;
            result.n = 2 * n;
            result.h = (b - a) / result.n;
            result.valueH = valueH;
            result.valueHalfH = valueHalfH;
            result.rungeCorrection = rungeCorrection;
            result.integralValue = integralValue;
            return result;
        }

        n *= 2;
    }
}

void printResult(const string& name, double eps, const Result& result)
{
    cout << name << endl;
    cout << "eps = " << eps << endl;
    cout << "n = " << result.n << endl;
    cout << "h = " << result.h << endl;
    cout << "I_h = " << result.valueH << endl;
    cout << "I_h/2 = " << result.valueHalfH << endl;
    cout << "Runge correction = " << result.rungeCorrection << endl;
    cout << "Integral = " << result.integralValue << endl;
    cout << endl;
}

int main()
{
    const double a = 0.0;
    const double b = 1.0;

    double eps;

    cout << "Enter eps: ";
    cin >> eps;

    if (eps <= 0.0)
    {
        cout << "Error: eps must be positive" << endl;
        return 1;
    }

    cout << fixed << setprecision(10);

    cout << "==============================" << endl;
    cout << "eps = " << eps << endl << endl;

    Result rect = solveByRunge(rectangleFormula, 3.0, 1, a, b, eps);
    Result trap = solveByRunge(trapezoidFormula, 3.0, -1, a, b, eps);
    Result simp = solveByRunge(simpsonFormula, 15.0, -1, a, b, eps);

    printResult("Rectangle formula", eps, rect);
    printResult("Trapezoid formula", eps, trap);
    printResult("Simpson formula", eps, simp);

    return 0;
}