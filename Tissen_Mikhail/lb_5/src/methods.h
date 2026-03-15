#ifndef METHODS_H
#define METHODS_H

#include <cmath>

// Внешние функции варианта задаются в main.cpp.
extern double F(double x);
extern double DF(double x);
extern double DDF(double x);

// Установка точности округления значений функции и производной.
// Если Delta = 0, округление не выполняется.
void SetRoundDelta(double Delta);

// Функция Round - округление числа X до ближайшего значения, кратного Delta.
double Round(double X, double Delta);

// Установка констант для критерия остановки метода Ньютона.
void SetNewtonConstants(double m1, double M2);

// Метод Ньютона для поиска корня уравнения F(x) = 0.
// X   - начальное приближение корня
// Eps - требуемая погрешность вычисления корня
// N   - число выполненных итераций
// Округление значений функции и производной задается через SetRoundDelta.
// Критерий остановки: |x_n - x_(n-1)| < sqrt(2*m1*Eps/M2)
double NEWTON(double X, double Eps, int &N);

#endif
