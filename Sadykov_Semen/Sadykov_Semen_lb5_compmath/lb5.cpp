#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Глобальная переменная для точности округления
double CurrentDelta = 0.0;

// Прототипы функций
double Round(double X, double Delta);
double F(double x);
double F1(double x);
double NEWTON(double X, double Eps, int &N);

// Реализация функций
double Round(double X, double Delta) {
    if (Delta <= 1E-9) return X;
    if (X > 0.0) return (Delta * (long((X / Delta) + 0.5)));
    else return (Delta * (long((X / Delta) - 0.5)));
}

double F(double x) {
    double raw = exp(x) - acos(sqrt(x));
    if (CurrentDelta > 1e-9) return Round(raw, CurrentDelta);
    return raw;
}

double F1(double x) {
    // Производная: e^x + 1 / (2 * sqrt(x - x^2))
    double raw = exp(x) + 1.0 / (2.0 * sqrt(x - x * x));
    if (CurrentDelta > 1e-9) return Round(raw, CurrentDelta);
    return raw;
}

double NEWTON(double X, double Eps, int &N) {
    double Y, Y1, DX;
    N = 0;
    do {
        Y = F(X);
        if (Y == 0.0) return X;
        Y1 = F1(X);
        if (fabs(Y1) < 1e-12) exit(1); // Защита от деления на 0
        DX = Y / Y1;
        X = X - DX;
        N++;
    } while (fabs(DX) > Eps && N < 100);
    return X;
}

int main() {
    // Начальное приближение по условию Фурье (левая граница)
    double x0 = 0.154405;

    // --- БЛОК 1: Зависимость N от Eps ---
    double EpsValues[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
    int numEps = 6;

    // Выводим только количество элементов для первого цикла Python
    printf("%d\n", numEps);
    
    CurrentDelta = 0.0; // Отключаем округление для первой части
    for (int i = 0; i < numEps; i++) {
        int iter;
        double res = NEWTON(x0, EpsValues[i], iter);
        // Формат: eps n x (разделенные пробелом)
        printf("%.10f %d %.10f\n", EpsValues[i], iter, res);
    }

    // --- БЛОК 2: Чувствительность к Delta ---
    double DeltaValues[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    int numDelta = 6;
    double fixedEps = 1e-6;

    // Выводим количество элементов для второго цикла Python
    printf("%d\n", numDelta);
    
    for (int i = 0; i < numDelta; i++) {
        CurrentDelta = DeltaValues[i]; // Включаем округление
        int iter;
        double res = NEWTON(x0, fixedEps, iter);
        // Формат: delta n x
        printf("%.10f %d %.10f\n", DeltaValues[i], iter, res);
    }

    return 0;
}