#ifndef METHODS_H
#define METHODS_H

// Внешняя функция F(x) - определяется в main.cpp для конкретного варианта
extern double F(double);

// Функция округления значения X до кратного Delta
// Используется для моделирования ошибок вычислений функции
double Round(double X, double Delta);

// Метод бисекции для поиска корня F(x) = 0
// Left, Right - границы начального интервала
// Eps - требуемая точность корня
// N - число итераций (возвращается по ссылке)
// Delta - точность округления значений функции (0 = без округления)
double BISECT(double Left, double Right, double Eps, int &N, double Delta = 0.0);

#endif