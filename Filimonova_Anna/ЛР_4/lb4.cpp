#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   

double CurrentDelta = 0.0;

double Round(double X, double Delta)
{
    if (Delta <= 1E-9)
    {
        fprintf(stderr, "Ошибка: Delta слишком мала\n");
        exit(1);
    }
    if (X > 0.0)
        return Delta * (long)((X / Delta) + 0.5);
    else
        return Delta * (long)((X / Delta) - 0.5);
}

double F(double x)
{
    double raw = exp(-x) - x * x * x;
    if (CurrentDelta > 0.0)
        return Round(raw, CurrentDelta);
    return raw;
}

double F1(double x)
{
    double raw = -exp(-x) - 3.0 * x * x;
    if (CurrentDelta > 0.0)
        return Round(raw, CurrentDelta);
    return raw;
}

double Phi(double x)
{
    double raw = exp(-x / 3.0);
    if (CurrentDelta > 0.0)
        return Round(raw, CurrentDelta);
    return raw;
}

double BISECT(double Left, double Right, double Eps, int &N)
{
    double E = fabs(Eps) * 2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left + Right) / 2.0;
    double Y;

    if (FLeft * FRight > 0.0) { puts("Неверное задание интервала\n"); exit(1); }
    if (Eps <= 0.0) { puts("Неверное задание точности\n"); exit(1); }

    N = 0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;

    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left);     
        Y = F(X);
        if (Y == 0.0) return X;
        if (Y * FLeft < 0.0)
            Right = X;
        else
        {
            Left = X;
            FLeft = Y;
        }
        N++;
    };
    return X;
}

double HORDA(double Left, double Right, double Eps, int &N)
{
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;

    if (FLeft * FRight > 0.0) { puts("Неверное задание интервала\n"); exit(1); }
    if (Eps <= 0.0) { puts("Неверное задание точности\n"); exit(1); }

    N = 0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;

    do
    {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        if (Y == 0.0) return X;
        if (Y * FLeft < 0.0)
        {
            Right = X;
            FRight = Y;
        }
        else
        {
            Left = X;
            FLeft = Y;
        }
        N++;
    }
    while (fabs(Y) >= Eps);
    return X;
}

double NEWTON(double X0, double Eps, int &N)
{
    double X = X0;
    double Y, Y1, DX;
    N = 0;
    do
    {
        Y = F(X);
        if (Y == 0.0) return X;

        Y1 = F1(X);
        if (Y1 == 0.0) { puts("Производная обратилась в ноль\n"); exit(1); }

        DX = Y / Y1;
        X = X - DX;
        N++;
    }
    while (fabs(DX) > Eps);
    return X;
}

double ITER(double X0, double Eps, int &N)
{
    if (Eps <= 0.0) { puts("Неверное задание точности\n"); exit(1); }
    double X1 = Phi(X0);
    double X2 = Phi(X1);
    N = 2;
    while ((X1 - X2) * (X1 - X2) > fabs((2 * X1 - X0 - X2) * Eps))
    {
        X0 = X1;
        X1 = X2;
        X2 = Phi(X1);
        N++;
    }
    return X2;
}

void interactiveMode()
{
    printf("e^{-x} - x^3 = 0\n");
    printf("Выберите метод:\n");
    printf("1 - Бисекция\n");
    printf("2 - Хорд\n");
    printf("3 - Ньютон\n");
    printf("4 - Простые итерации\n");
    printf("Ваш выбор: ");
    int choice;
    scanf("%d", &choice);

    double left = 0.5, right = 1.0;
    double x0 = 1.0;
    double eps, delta;
    int studyType;
    printf("Тип исследования:\n");
    printf("1 - Зависимость от точности Eps\n");
    printf("2 - Чувствительность к ошибкам округления Delta\n");
    scanf("%d", &studyType);

    if (studyType == 1)
    {
        printf("Точность Eps: ");
        scanf("%lf", &eps);
        CurrentDelta = 0.0;
        int iter;
        double root;
        switch (choice)
        {
        case 1: root = BISECT(left, right, eps, iter); break;
        case 2: root = HORDA(left, right, eps, iter); break;
        case 3: root = NEWTON(x0, eps, iter); break;
        case 4: root = ITER(x0, eps, iter); break;
        default: printf("Неверный выбор\n"); return;
        }
        printf("Корень = %.10f, число итераций = %d\n", root, iter);
    }
    else if (studyType == 2)
    {
        printf("Фиксированная точность Eps: ");
        scanf("%lf", &eps);
        printf("Величина ошибки округления Delta: ");
        scanf("%lf", &delta);
        CurrentDelta = delta;
        int iter;
        double root;
        switch (choice)
        {
        case 1: root = BISECT(left, right, eps, iter); break;
        case 2: root = HORDA(left, right, eps, iter); break;
        case 3: root = NEWTON(x0, eps, iter); break;
        case 4: root = ITER(x0, eps, iter); break;
        default: printf("Неверный выбор\n"); return;
        }
        printf("Корень с округлением = %.10f, число итераций = %d\n", root, iter);
    }
    else
    {
        printf("Неверный тип исследования\n");
    }
}

void printHelp(const char* progname)
{
    printf("format:\n");
    printf("%s method type left right x0 fixed_eps [params...]\n", progname);
    printf("\n");
    printf("  method      - bisect, horda, newton, iter\n");
    printf("  type        - eps, delta\n");
    printf("  left, right - границы интервала (bisect, horda)\n");
    printf("  x0          - начальное приближение (newton, iter)\n");
    printf("  fixed_eps   - фиксированная точность для delta\n");
    printf("  параметры   - список значений eps или delta\n");
    printf("\n");
}

int main(int argc, char *argv[])
{

    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        printHelp(argv[0]);
        return 0;
    }

    if (argc == 1)
    {
        interactiveMode();
        return 0;
    }

    if (argc < 8)
    {
        fprintf(stderr, "Ошибка: недостаточно аргументов. Используйте -h для справки.\n");
        return 1;
    }

    char* method_str = argv[1];
    char* type_str = argv[2];
    double left = atof(argv[3]);
    double right = atof(argv[4]);
    double x0 = atof(argv[5]);
    double fixed_eps = atof(argv[6]);

    int method_id = -1;
    if (strcmp(method_str, "bisect") == 0) method_id = 1;
    else if (strcmp(method_str, "horda") == 0) method_id = 2;
    else if (strcmp(method_str, "newton") == 0) method_id = 3;
    else if (strcmp(method_str, "iter") == 0) method_id = 4;
    else
    {
        fprintf(stderr, "Неизвестный метод: %s\n", method_str);
        return 1;
    }

    int type_id = -1;
    if (strcmp(type_str, "eps") == 0) type_id = 1;
    else if (strcmp(type_str, "delta") == 0) type_id = 2;
    else
    {
        fprintf(stderr, "Неизвестный тип: %s\n", type_str);
        return 1;
    }

    int nParams = argc - 7;
    double* params = new double[nParams];
    for (int i = 0; i < nParams; i++)
    {
        params[i] = atof(argv[7 + i]);
    }

    for (int i = 0; i < nParams; i++)
    {
        double param = params[i];
        int iter;
        double root;

        if (type_id == 2)
            CurrentDelta = param;
        else
            CurrentDelta = 0.0;

        if (method_id == 1)
            root = BISECT(left, right, (type_id == 1) ? param : fixed_eps, iter);
        else if (method_id == 2)
            root = HORDA(left, right, (type_id == 1) ? param : fixed_eps, iter);
        else if (method_id == 3)
            root = NEWTON(x0, (type_id == 1) ? param : fixed_eps, iter);
        else
            root = ITER(x0, (type_id == 1) ? param : fixed_eps, iter);

        printf("%.15e %d %.15f\n", param, iter, root);
    }

    delete[] params;
    return 0;
}