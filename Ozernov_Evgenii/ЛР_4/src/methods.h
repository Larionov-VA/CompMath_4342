#pragma once
#include <functional>

struct HordaResult {
    double x;        // приближение корня
    double a, b;     // текущий отрезок (сохраняем "скобку")
    int iters;       // число итераций
    bool ok;         // успешно ли сошлось
};

double Round(double value, double Delta);

// Метод хорд (regula falsi):
// f      - функция
// Left,Right - начальный отрезок, где f(Left)*f(Right)<0
// Eps    - точность
// nmax   - лимит итераций
// Delta  - округление значений f(...) через Round (0 => без округления)
HordaResult HORDA(const std::function<double(double)>& f,
                  double Left, double Right,
                  double Eps,
                  int nmax = 1000000,
                  double Delta = 0.0);