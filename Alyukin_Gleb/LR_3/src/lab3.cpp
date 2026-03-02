#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"

double delta; 

double F(double x)
{
    extern double delta;
    double s;

    s = x * x - 5.0 * sin(x); // Значение функции

    s = Round(s, delta); // Округление с точностью delta

    return s;
}

int main(){
    double a, b, eps, x;
    int k;
    float a1, b1, eps1, delta1;

    printf("Лабораторная работа №3. Вариант 1\n");
    printf("Функция: f(x) = x^2 - 5*sin(x)\n\n");

    a = 2;
    b = 3;
    eps = 0.000001; //Точность вычисления корня
    delta = 0.0000001; // Точность задания функции

    x = BISECT(a, b, eps, k);

    printf("\nРезультаты:\n");
    printf("Корень x = %f\n", x);
    printf("Число итераций N = %d\n", k);
    return 0;
}
