#include <iostream>
#include "methods.hpp"
#include <math.h>
#include <fstream>

double DELTA = 1E-16;

double F(double x)
{
    return Round((1.0+std::cos(x))/(3.0-std::sin(x)) - x, DELTA);
}


void study_operation_count(){
    std::ofstream out;
    out.open("./operation_metrix.txt");

    double p = 0.1;

    if (out.is_open())
    {
        for (int i = 0; i < 20; i++)
        {
            int N = 0;
            double root = BISECT(0.0, 1.0, p, N);

            out << p << ' ' << N << ' ' << root << std::endl;
            p /= 10.0;
        }
    }

    out.close();
}

void study_root_value(){
    
    std::ofstream out;
    out.open("./root_metrix.txt");

    DELTA = 0.1;

    if (out.is_open())
    {
        for (int i = 0; i < 10; i++){
            int N = 0;

            double root = BISECT(0.0, 1.0, 0.01, N);

            out << root << ' ' << N << ' ' << DELTA << std::endl;

            DELTA /= 10;
        }
    }
    out.close();
}

int main(){

    study_operation_count();

    study_root_value();

    return 0;
}