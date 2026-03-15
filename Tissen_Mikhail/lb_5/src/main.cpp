#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"

// ============================================================
// ВАРИАНТ 18: f(x) = e^(1/x^2) - ln(x)
// Область определения: x > 0
// Интервал изоляции: [3; 4]
// ============================================================
double F(double x)
{
    if (x <= 0.0) {
        printf("Ошибка: x должен быть > 0\n");
        exit(1);
    }
    return exp(1.0 / (x * x)) - log(x);
}

double DF(double x)
{
    if (x <= 0.0) {
        printf("Ошибка: x должен быть > 0\n");
        exit(1);
    }
    return -2.0 * exp(1.0 / (x * x)) / (x * x * x) - 1.0 / x;
}

double DDF(double x)
{
    if (x <= 0.0) {
        printf("Ошибка: x должен быть > 0\n");
        exit(1);
    }
    return (pow(x, 4.0) + 6.0 * x * x * exp(1.0 / (x * x)) + 4.0 * exp(1.0 / (x * x))) / pow(x, 6.0);
}

// ============================================================
// Головная программа
// ============================================================
int main()
{
    double Left = 3.0;
    double Right = 4.0;
    // Начальное приближение выбираем из условия f(x0) * f''(x0) > 0.
    double X0 = 3.0;
    double RefRoot, Root;
    double Eps, Delta;
    double m1, M2;
    int N;

    printf("========================================\n");
    printf("  Laboratornaya rabota #5: Metod Newtona\n");
    printf("  Variant 18: f(x) = e^(1/x^2) - ln(x)\n");
    printf("========================================\n\n");

    printf("Interval otdeleniya kornya: [%.2f; %.2f]\n", Left, Right);
    printf("F(Left)   = %.10f\n", F(Left));
    printf("F(Right)  = %.10f\n", F(Right));
    printf("DF(Left)  = %.10f\n", DF(Left));
    printf("DF(Right) = %.10f\n", DF(Right));
    printf("DDF(Left) = %.10f\n", DDF(Left));
    printf("DDF(Right)= %.10f\n\n", DDF(Right));

    // Оценки min |f'(x)| и max |f''(x)| на отрезке [3, 4].
    m1 = fabs(DF(Right));
    M2 = fabs(DDF(Left));
    SetNewtonConstants(m1, M2);
    printf("m1 = min |f'(x)| = %.10f\n", m1);
    printf("M2 = max |f''(x)| = %.10f\n", M2);
    printf("Nachalnoe priblizhenie x0 = %.10f\n\n", X0);

    // ============================================================
    // Нахождение опорного корня с высокой точностью
    // ============================================================
    SetRoundDelta(0.0);
    RefRoot = NEWTON(X0, 1e-12, N);
    printf("Opornyy koren (Eps=1e-12, Delta=0): %.15f (iteraciy: %d)\n\n", RefRoot, N);

    // ============================================================
    // Исследование скорости сходимости
    // ============================================================
    printf("Zavisimost chisla iteraciy ot Eps:\n");
    printf("Eps\t\tIteracii\tKoren\n");
    printf("-----------------------------------------------\n");

    double EpsValues[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12};
    int nEps = sizeof(EpsValues) / sizeof(EpsValues[0]);

    FILE* feps = fopen("newton_iterations_vs_eps.csv", "w");
    fprintf(feps, "Eps,Iterations,Root\n");

    for (int i = 0; i < nEps; i++) {
        Eps = EpsValues[i];
        SetRoundDelta(0.0);
        Root = NEWTON(X0, Eps, N);
        printf("%.1e\t%d\t\t%.15f\n", Eps, N, Root);
        fprintf(feps, "%e,%d,%.15f\n", Eps, N, Root);
    }
    fclose(feps);
    printf("\nDannye sohraneny v newton_iterations_vs_eps.csv\n");

    // ============================================================
    // Исследование чувствительности к ошибкам округления
    // ============================================================
    printf("\nChuvstvitelnost metoda (Eps=1e-8):\n");
    printf("Delta\t\tIteracii\tKoren\t\t\tError\n");
    printf("--------------------------------------------------------------------------\n");

    double DeltaValues[] = {0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
    int nDelta = sizeof(DeltaValues) / sizeof(DeltaValues[0]);

    FILE* fdel = fopen("newton_sensitivity_vs_delta.csv", "w");
    fprintf(fdel, "Delta,Iterations,Root,Error\n");

    for (int i = 0; i < nDelta; i++) {
        Delta = DeltaValues[i];
        SetRoundDelta(Delta);
        Root = NEWTON(X0, 1e-8, N);
        double Err = fabs(Root - RefRoot);
        printf("%.1e\t%d\t\t%.15f\t%.15f\n", Delta, N, Root, Err);
        fprintf(fdel, "%e,%d,%.15f,%.15f\n", Delta, N, Root, Err);
    }
    fclose(fdel);
    printf("\nDannye sohraneny v newton_sensitivity_vs_delta.csv\n");

    return 0;
}
