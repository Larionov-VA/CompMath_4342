#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"

double delta;

double F(double x){
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

    printf("Лабораторная работа №4. Вариант 1\n");
    printf("Функция: f(x) = x^2 - 5*sin(x)\n\n");

    printf("Введите границы интервала a и b: "); // 2 3
    scanf("%f %f", &a1, &b1);
    a = a1;
    b = b1;

    printf("Введите точность вычисления корня: "); // 0.0001
    scanf("%f", &eps1);
    eps = eps1; //Точность вычисления корня

    printf("Введите точность задания функции: "); // 0.000001
    scanf("%f", &delta1);
    delta = delta1; // Точность задания функции

    x = HORDA(a, b, eps, k);

    printf("\nРезультаты:\n");
    printf("x = %f\n", x);
    printf("Число итераций N = %d\n", k);
    return 0;
}
