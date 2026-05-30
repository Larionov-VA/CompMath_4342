#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"


// f(x) = sin(x^2) - 6x + 1
// Интервал изоляции: [0; 0.5]

double F(double x)
{
    return sin(x*x) - 6*x + 1;
}

double DF(double x)
{
    return 2*x*cos(x*x) - 6;
}

double DDF(double x)
{
    return 2*cos(x) - 2*x*sin(x);
}


int main()
{
    double Left = 0.0;
    double Right = 0.5;
    // Начальное приближение выбираем из условия f(x0) * f''(x0) > 0.
    double X0 = 0.0;
    double RefRoot, Root;
    double Eps, Delta;
    double m1, M2;
    int N;

    printf("f(x) = sin(x^2) - 6x + 1\n");

    printf("Interval otdeleniya kornya: [%.2f; %.2f]\n", Left, Right);
    printf("F(Left)   = %.10f\n", F(Left));
    printf("F(Right)  = %.10f\n", F(Right));
    printf("DF(Left)  = %.10f\n", DF(Left));
    printf("DF(Right) = %.10f\n", DF(Right));
    printf("DDF(Left) = %.10f\n", DDF(Left));
    printf("DDF(Right)= %.10f\n\n", DDF(Right));

    // Оценки min |f'(x)| и max |f''(x)| на отрезке [0, 0.5].
    m1 = fabs(DF(Right));
    M2 = fabs(DDF(Left));
    SetNewtonConstants(m1, M2);
    printf("m1 = min |f'(x)| = %.10f\n", m1);
    printf("M2 = max |f''(x)| = %.10f\n", M2);
    printf("x0 = %.10f\n\n", X0);

    // Нахождение опорного корня с высокой точностью
    SetRoundDelta(0.0);
    RefRoot = NEWTON(X0, 1e-12, N);
    printf("(Eps=1e-12, Delta=0): %.15f (Iterations: %d)\n\n", RefRoot, N);

    // Исследование скорости сходимости
    printf("Eps\t\tIterations\tRoot\n");
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

    // Исследование чувствительности к ошибкам округления
    printf("\n(Eps=1e-8):\n");
    printf("Delta\t\tIterations\tRoot\t\t\tError\n");
    
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

    return 0;
}