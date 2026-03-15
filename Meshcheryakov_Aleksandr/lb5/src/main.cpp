#include <iostream>
#include <fstream>
#include <cmath>
#include "../../methods.hpp"

// Метод Ньютона

double DELTA = 1E-16;

double F(double x)
{
    return Round((1.0 + std::cos(x)) / (3.0 - std::sin(x)) - x, DELTA);
}

double F1(double x){
    return Round(((std::cos(x)*(std::cos(x) + 1.0) - (3.0 - std::sin(x))*std::sin(x)) / std::pow((3.0 - std::sin(x)), 2)) - 1, DELTA);
}

int main()
{
    std::ofstream main_out;
    std::ofstream conditionality_out;

    main_out.open("main_metrix.txt");
    int N;

    if (main_out.is_open())
    {
        for (double eps = 0.1; eps >= 1E-15; eps /= 10)
        {
            double result = NEWTON((1.0 - 0.0)/2.0, eps, N);
            main_out << result << ' ' << eps << ' ' << N << std::endl;
        }
    }
    main_out.close();

    DELTA = 0.1;
    conditionality_out.open("conditionality_metrix.txt");

    if (conditionality_out.is_open())
    {
        for (int i = 0; i < 10; i++)
        {
            double result = NEWTON((1.0 - 0.0) / 2.0, 1E-15, N);
            conditionality_out << result << ' ' << DELTA << ' ' << N << std::endl;
            DELTA /= 10;
        }
    }
    conditionality_out.close();

    return 0;
}
