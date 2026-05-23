#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Глобальная переменная для точности округления
double CurrentDelta = 0.0;

// Прототипы функций
double Round(double X, double Delta);
double PHI(double x); // Переименовали для ясности, так как это x = phi(x)
double ITER(double X0, double Eps, int &N);

// Функция округления
double Round(double X, double Delta) {
    if (Delta <= 1E-9) return X;
    if (X > 0.0) return (Delta * (long((X / Delta) + 0.5)));
    else return (Delta * (long((X / Delta) - 0.5)));
}

// Итерационная функция phi(x) = cos^2(e^x)
double PHI(double x) {
    double raw = pow(cos(exp(x)), 2);
    if (CurrentDelta > 1e-9) return Round(raw, CurrentDelta);
    return raw;
}

// Метод простых итераций
double ITER(double X0, double Eps, int &N) {
    if (Eps <= 0.0) exit(1);
    
    double X1 = PHI(X0);
    double X2 = PHI(X1);
    N = 2;
    
    // Условие выхода
    while (fabs(X1 - X2) > Eps) {
        X0 = X1;
        X1 = X2;
        X2 = PHI(X1);
        N++;
        if (N > 1000) break; // Защита от бесконечного цикла, если условие сжатия не выполнено
    }
    return X2;
}

int main() {
    double x0 = 0.1; // Начальное приближение

    // --- БЛОК 1: Зависимость N от Eps ---
    double EpsValues[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
    int numEps = 6;

    printf("%d\n", numEps); // Кол-во строк для первого цикла Python
    CurrentDelta = 0.0;
    for (int i = 0; i < numEps; i++) {
        int iter;
        double res = ITER(x0, EpsValues[i], iter);
        // Формат: eps n x
        printf("%.10f %d %.10f\n", EpsValues[i], iter, res);
    }

    // --- БЛОК 2: Чувствительность к Delta ---
    double DeltaValues[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    int numDelta = 6;
    double fixedEps = 1e-6;

    printf("%d\n", numDelta); // Кол-во строк для второго цикла Python
    for (int i = 0; i < numDelta; i++) {
        CurrentDelta = DeltaValues[i];
        int iter;
        double res = ITER(x0, fixedEps, iter);
        // Формат: delta n x
        printf("%.10f %d %.10f\n", DeltaValues[i], iter, res);
    }

    return 0;
}