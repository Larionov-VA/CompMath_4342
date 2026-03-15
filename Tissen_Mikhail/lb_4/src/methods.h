#pragma once

#include <functional>

struct HordaResult {
    double x;      // найденное приближение корня
    double a, b;   // текущий отрезок, сохраняющий смену знака
    int iters;     // число выполненных итераций
    bool ok;       // достигнут ли критерий остановки
};

// Функция Round - округление числа до кратного delta.
// Используется для моделирования ошибок вычисления значений функции.
double Round(double value, double delta);

// Метод хорд для поиска корня уравнения f(x) = 0.
// f            - функция одной переменной
// left, right  - границы начального отрезка изоляции
// eps          - требуемая точность корня
// nmax         - предельное число итераций
// delta        - шаг округления значений f(x) (0 = без округления)
HordaResult HORDA(const std::function<double(double)>& f,
                  double left,
                  double right,
                  double eps,
                  int nmax = 1000000,
                  double delta = 0.0);
