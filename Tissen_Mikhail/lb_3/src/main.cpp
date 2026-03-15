#include <stdio.h>
#include <math.h>
#include "methods.h"

// ============================================================
// ВАРИАНТ 18: f(x) = e^(1/x²) - ln(x)
// Область определения: x > 0
// Интервал изоляции: [3; 4] (f(3) > 0, f(4) < 0)
// ============================================================
double F(double x)
{
    if (x <= 0.0) {
        printf("Ошибка: x должен быть > 0\n");
        exit(1);
    }
    return exp(1.0 / (x * x)) - log(x);
}

// ============================================================
// Головная программа
// ============================================================
int main()
{
    // Интервал изоляции корня для Варианта 18: [3; 4]
    double Left = 3.0;
    double Right = 4.0;
    double Eps, Root, Delta, RefRoot;
    int N;

    printf("========================================\n");
    printf("  Laboratornaya rabota #3: Metod bisekcii\n");
    printf("  Variant 18: f(x) = e^(1/x^2) - ln(x)\n");
    printf("========================================\n\n");

    // Вывод интервала и проверка условия теоремы Коши
    printf("Interval isolyacii: [%.2f; %.2f]\n", Left, Right);
    printf("F(Left) = %.6f\n", F(Left));
    printf("F(Right) = %.6f\n", F(Right));
    printf("F(Left)*F(Right) = %.6f %s\n\n", 
           F(Left)*F(Right), 
           F(Left)*F(Right) < 0 ? "(uslovie vypolneno)" : "(OSHIBKA!)");

    // Нахождение опорного корня с высокой точностью (для сравнения)
    // Delta = 0 означает отсутствие округления (точное вычисление)
    RefRoot = BISECT(Left, Right, 1e-10, N, 0.0);
    printf("Opornyy koren (Eps=1e-10, Delta=0): %.10f (iteraciy: %d)\n\n", RefRoot, N);

    // ============================================================
    // Исследование зависимости числа итераций от точности Eps
    // Диапазон: от 0.1 до 0.000001 (6 значений)
    // ============================================================
    printf("Zavisimost iteraciy ot Eps:\n");
    printf("Eps\t\tIteracii\tKoren\n");
    printf("----------------------------------------\n");

    // Включаем Eps = 0.1 в массив значений
    double EpsValues[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
    int nEps = sizeof(EpsValues) / sizeof(EpsValues[0]);

    // Сохранение данных для построения графика
    FILE* feps = fopen("iterations_vs_eps.csv", "w");
    fprintf(feps, "Eps,Iterations,Root\n");

    for (int i = 0; i < nEps; i++) {
        Eps = EpsValues[i];
        // Delta = 0 - без округления для чистоты эксперимента
        Root = BISECT(Left, Right, Eps, N, 0.0);
        printf("%.1e\t%d\t\t%.10f\n", Eps, N, Root);
        fprintf(feps, "%e,%d,%.10f\n", Eps, N, Root);
    }
    fclose(feps);
    printf("\nDannye sohraneny: iterations_vs_eps.csv\n");

    // ============================================================
    // Исследование чувствительности к ошибкам (параметр Delta)
    // Фиксированная точность Eps, меняем точность округления функции
    // ============================================================
    printf("\nChuvstvitelnost metoda (Eps=1e-6):\n");
    printf("Delta\t\tKoren\t\tPogreshnost\n");
    printf("----------------------------------------\n");

    // Включаем Delta = 0 для сравнения
    double DeltaValues[] = {0.0, 0.1, 0.01, 0.001, 0.0001, 0.00001};
    int nDelta = sizeof(DeltaValues) / sizeof(DeltaValues[0]);

    FILE* fdel = fopen("sensitivity_vs_delta.csv", "w");
    fprintf(fdel, "Delta,Root,Error\n");

    for (int i = 0; i < nDelta; i++) {
        Delta = DeltaValues[i];
        // Округление применяется ВНУТРИ BISECT к значениям функции
        Root = BISECT(Left, Right, 1e-6, N, Delta);
        double Err = fabs(Root - RefRoot);
        printf("%.1e\t%.10f\t%.10f\n", Delta, Root, Err);
        fprintf(fdel, "%e,%.10f,%.10f\n", Delta, Root, Err);
    }
    fclose(fdel);
    printf("\nDannye sohraneny: sensitivity_vs_delta.csv\n");

    return 0;
}