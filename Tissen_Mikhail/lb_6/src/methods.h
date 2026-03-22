#pragma once

#include <functional>

struct IterResult {
    double x;         // последнее приближение
    double lastDiff;  // |x_(n+1) - x_n| на последнем шаге
    int iters;        // число выполненных итераций
    bool ok;          // true, если достигнут критерий остановки
};

// Округление значения до ближайшего числа, кратного delta.
// Если delta <= 0, значение возвращается без изменений.
double Round(double value, double delta);

// Оценка q = max |phi'(x)| на [left, right] равномерным перебором точек.
double EstimateQ(const std::function<double(double)>& dphi,
                 double left,
                 double right,
                 int samples = 10000);

// Метод простых итераций для уравнения x = phi(x).
// phi   - итерационная функция
// x0    - начальное приближение
// eps   - требуемая точность по условию |x_(n+1) - x_n|
// nmax  - максимально допустимое число итераций
IterResult ITER(const std::function<double(double)>& phi,
                double x0,
                double eps,
                int nmax = 1000);
