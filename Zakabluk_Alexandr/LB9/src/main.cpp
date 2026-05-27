#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#include <vector>
#include <string>
#include <windows.h> // Для кодировки консоли

/*
 Подынтегральная функция.
 По условию варианта 5:
 f(x) = exp(cos(x))
 На отрезке [0; 1]
*/
double f(double x) {
    return std::exp(std::cos(x));
}

/*
 Составная формула прямоугольников (средних).
 I ≈ h * сумма f(x_i + h/2)
*/
double rectangleFormula(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + (i + 0.5) * h; // Середина отрезка
        sum += f(x);
    }
    return h * sum;
}

/*
 Составная формула трапеций.
 I ≈ h * ( (f(a) + f(b)) / 2 + сумма f(x_i) )
*/
double trapezoidFormula(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }
    return h * sum;
}

/*
 Составная формула Симпсона.
 I ≈ h/3 * (f(x0) + 4*f(x1) + 2*f(x2) + ... + f(xn))
 Внимание: n должно быть четным!
*/
double simpsonFormula(double a, double b, int n) {
    if (n % 2 != 0) ++n; // Гарантируем четность
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 1) {
            sum += 4.0 * f(x); // Нечетные индексы
        } else {
            sum += 2.0 * f(x); // Четные индексы
        }
    }
    return (h / 3.0) * sum;
}

/*
 Структура для хранения результата работы метода.
*/
struct Result {
    std::string methodName;
    double eps;
    int n;
    double h;
    double integral;
    double rungeIntegral;
    double error;
    int iterations;
};

/*
 Универсальная функция для применения правила Рунге.
*/
Result rungeRule(
    const std::string& methodName,
    const std::function<double(double, double, int)>& formula,
    int order, // Порядок точности (2 для прямоуг/трапеций, 4 для Симпсона)
    double a,
    double b,
    double eps
) {
    int n = 2; // Начинаем с n=2
    double oldIntegral = formula(a, b, n);
    int iterations = 0;

    while (true) {
        n *= 2; // Удваиваем сетку
        ++iterations;
        
        double newIntegral = formula(a, b, n);
        
        // Знаменатель в правиле Рунге: 2^p - 1
        double denominator = std::pow(2.0, order) - 1.0;
        
        // Оценка погрешности
        double error = std::abs(newIntegral - oldIntegral) / denominator;
        
        // Уточненное значение по Рунге (с учетом знака ошибки)
        double signedError = (newIntegral - oldIntegral) / denominator;
        double rungeIntegral = newIntegral + signedError;

        if (error <= eps) {
            Result result;
            result.methodName = methodName;
            result.eps = eps;
            result.n = n;
            result.h = (b - a) / n;
            result.integral = newIntegral;
            result.rungeIntegral = rungeIntegral;
            result.error = error;
            result.iterations = iterations;
            return result;
        }
        oldIntegral = newIntegral;
    }
}

// Функция вывода
void printResult(const Result& result) {
    std::cout << std::left << std::setw(20) << result.methodName
              << "| eps: " << std::setw(6) << result.eps
              << "| n: " << std::setw(3) << result.n
              << "| h: " << std::setw(8) << result.h
              << "| I = " << std::setw(12) << result.integral
              << "| I(Рунге) = " << std::setw(12) << result.rungeIntegral
              << "| Err: " << std::setw(12) << result.error << "\n";
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    
    std::cout << std::fixed << std::setprecision(10);

    // Параметры интегрирования
    double a = 0.0;
    double b = 1.0;
    std::vector<double> epsValues = {0.01, 0.001, 0.000001};

    std::cout << "Лабораторная работа №9\n";
    std::cout << "Интеграл: exp(cos(x)) dx на отрезке [0; 1]\n";
    std::cout << std::string(110, '-') << "\n\n";

    for (double eps : epsValues) {
        std::cout << ">>> Требуемая точность eps = " << eps << "\n";
        
        Result rRect = rungeRule("Прямоугольники", rectangleFormula, 2, a, b, eps);
        Result rTrap = rungeRule("Трапеции", trapezoidFormula, 2, a, b, eps);
        Result rSimp = rungeRule("Симпсон", simpsonFormula, 4, a, b, eps);
        
        printResult(rRect);
        printResult(rTrap);
        printResult(rSimp);
        std::cout << std::string(110, '-') << "\n";
    }

    return 0;
}