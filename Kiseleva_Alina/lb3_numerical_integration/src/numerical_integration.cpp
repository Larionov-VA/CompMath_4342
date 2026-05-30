#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

double f(double x){
    return std::cosh(x * x);
}

double rectangle(std::vector<double>& arr_x, double h, int n){
    double sum = 0.0;

    for (int i = 0; i < n; ++i){
        sum += f(arr_x[i] + h/2);
    }

    return h * sum;
}

double trapezoid(std::vector<double>& arr_x, double h, int n){
    double sum = 0.0;

    for (int i = 0; i < n; ++i){
        sum += f(arr_x[i]) + f(arr_x[i+1]);
    }

    return h/2 * sum;
}

double simpson(std::vector<double>& arr_x, double h, int m){
    double sum = 0.0;

    for (int i = 0; i < m; ++i){
        sum += f(arr_x[2*i]) + 4*f(arr_x[2*i+1]) + f(arr_x[2*i+2]);
    }

    return h/3 * sum;
}

std::vector<double> make_x_net(double h, int n, double a){
    std::vector<double> arr_x(n + 1);
    for (int i = 0; i <= n; ++i){
        arr_x[i] = a + i*h;
    }
    return arr_x;
}

void runge(int method, int n, double a, double b, double eps){
    std::vector<std::pair<int, double>> ans;
    double runge_error = 1.0;
    double I = 0.0;
    int cnt_incr_n = 0;
    std::string method_name = "";

    switch(method){
        case 1: {
            method_name = "Rectangle";
            break;
        }
        case 2: {
            method_name = "Trapezoid";
            break;
        }
        case 3: {
            method_name = "Simpson";
            break;
        }
        default: {
                std::cout << "Неизвестный метод\n";
                return;
        }
    }

    std::cout << "\nМетод " << method_name << ":\n";

    while (runge_error >= eps){
        int m = n/2;
        double h = (b - a) / n;

        std::vector<double> arr_x = make_x_net(h, n, a);
        std::vector<double> arr_x_2 = make_x_net(h/2, 2*n+1, a);

        switch(method){
            case 1: {
                double I_rectangle_h = rectangle(arr_x, h, n);
                double I_rectangle_h_2 = rectangle(arr_x_2, h/2, 2*n);
                runge_error = fabs(I_rectangle_h_2 - I_rectangle_h) / 3.0;
                I = I_rectangle_h_2 + runge_error;
                break;
            }
            case 2: {
                double I_trapezoid_h = trapezoid(arr_x, h, n);
                double I_trapezoid_h_2 = trapezoid(arr_x_2, h/2, 2*n);
                runge_error = fabs(I_trapezoid_h_2 - I_trapezoid_h) / 3.0;
                I = I_trapezoid_h_2 - runge_error;
                break;
            }
            case 3: {
                double I_simpson_h = simpson(arr_x, h, m);
                double I_simpson_h_2 = simpson(arr_x_2, h/2, n);
                runge_error = fabs(I_simpson_h_2 - I_simpson_h) / 15.0;
                I = I_simpson_h_2 - runge_error;
                break;
            }
        }

        std::cout << "n = " << n << ", I = " << I << ", error = " << runge_error << '\n';
        n *= 2;
        cnt_incr_n++;
    }

    std::cout << "Количество удвоений: " << cnt_incr_n-1 << '\n';
}

int main(){
    std::cout << std::fixed << std::setprecision(12);

    double eps;
    std::cout << "Введите значение прогрешности epsilon:\n";
    std::cin >> eps;

    double a = 0.0;
    double b = 1.0;

    runge(1, 2, a, b, eps);
    runge(2, 2, a, b, eps);
    runge(3, 2, a, b, eps);

    return 0;
}