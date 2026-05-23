#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"
#define __NEWTON

double delta;

double F(double x){
    double res = x * x - 5.0 * sin(x);
    return Round(res, delta);
}

double derivativeF(double x) {
    double res = 2.0 * x - 5.0 * cos(x);
    return Round(res, delta);
}

int main(){
     double x0, eps, x;
    int k;
    float x0_in, eps1, delta1;

    printf("Лабораторная работа №5. Вариант 1\n");
    printf("Функция: f(x) = x^2 - 5*sin(x)\n");

    printf("Введите начальное приближение x0: "); // 3.0
    scanf("%f", &x0_in);
    x0 = 3.0;

    printf("Введите точность вычисления корня (Eps): "); // 0.0001
    scanf("%f", &eps1);
    eps = eps1;

    printf("Введите точность задания функции (Delta): "); // 0.00000001
    scanf("%f", &delta1);
    delta = delta1;

    x = NEWTON(x0, eps, k);

    printf("Результаты:\n");
    printf("x = %7.7f\n", x);
    printf("Число итераций N = %d\n", k);
    return 0;
}
