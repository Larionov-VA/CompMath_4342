#include <cstddef>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

double Eps = 0.000001;

double F(double x){
    return std::log(x) * (1.0/(x+1));
}

double Rectangle(double a, double b, int n){
    
    double h = (b - a)/n;
    double Sum = 0.0;

    for (int i = 0; i < n; i++){
        Sum += F((a + i*h) + (h/2.0));
    }

    return h * Sum;
}

double Trapezoids(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (F(a) + F(b)) / 2.0;
    for (int i = 1; i < n; i++) {
        sum += F(a + i * h);
    }
    return h * sum;
}


double Simpson(double a, double b, int n) {
    if (n % 2 != 0) n++; 
    double h = (b - a) / n;
    double sum = F(a) + F(b);
    for (int i = 1; i < n; i++) {
        if (i % 2 == 0) sum += 2.0 * F(a + i * h);
        else sum += 4.0 * F(a + i * h);
    }
    return (h / 3.0) * sum;
}

void reserch(int methods){ 
    int n = 2;
    double a = 1.0, b = 2.0;
    double I_h, I_h2, error;

    switch(methods){
        
        case(0):
            std::cout << "Rectangle" << std::endl;
            do {
                I_h = Rectangle(a, b, n);
                I_h2 = Rectangle(a, b, 2 * n);
                error = std::abs(I_h2 - I_h) / 3.0;
    
                if (error >= Eps) n *= 2;
            } while (error >= Eps);

            std::cout << "Result: " << I_h2 << " n: " << 2*n << std::endl;
            break;

        case(1):

            std::cout << "\nTrapezoid" << std::endl;
            do {
                I_h = Trapezoids(a, b, n);
                I_h2 = Trapezoids(a, b, 2 * n);
                error = std::abs(I_h2 - I_h) / 3.0;
    
                if (error >= Eps) n *= 2;
            } while (error >= Eps);

            std::cout << "Result: " << I_h2 << " n: " << 2*n << std::endl;
            break;

        case(2):
            std::cout << "\nSimpson" << std::endl; 
            do {
                I_h = Simpson(a, b, n);
                I_h2 = Simpson(a, b, 2 * n);
                error = std::abs(I_h2 - I_h) / 15.0;
    
                if (error >= Eps) n *= 2;
            } while (error >= Eps);
            

            int ep = std::to_string(Eps).size() - 2;

            std::cout << std::fixed << std::setprecision(ep);
            std::cout << "Result: " << I_h2 << " n: " << 2*n << std::endl;
            break;

    }

}

int main(){

    double e;
    std::cout << "Введите желаемую точность (eps): ";
    std::cin >> e; std::cout << std::endl;
    
    Eps = e;


    reserch(0);
    reserch(1);
    reserch(2);

    return 0;
}
