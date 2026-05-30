#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

double f(double x){
    return tan(x) - 1/x;
}

std::vector<double> solve_tridiagonal_matrix(std::vector<std::vector<double>> matrix, int n){
    int m = n - 2;

    std::vector<double> a(m);
    std::vector<double> b(m);
    std::vector<double> S(n, 0.0);

    a[0] = -matrix[0][1] / matrix[0][0];
    b[0] = matrix[0][m] / matrix[0][0];

    for (int i = 1; i < m-1; ++i){
        double A = matrix[i][i-1];
        double B = matrix[i][i];
        double C = matrix[i][i+1];
        double D = matrix[i][m];

        double denom = A*a[i-1] + B;
        
        a[i] = -C/denom;
        b[i] = (D - A*b[i-1])/denom;
    }

    S[m] = b[m-1];

    for (int i = m-2; i >= 0; --i)
        S[i+1] = a[i]*S[i+2] + b[i];

    S[0] = S[n-1] = 0.0;
    
    return S;
}

void spline_interpolation(std::vector<double> h, std::vector<double> y, int n, std::ofstream& f){
    int m = n - 2;
    std::vector<std::vector<double>> system_matrix(m, std::vector<double>(m + 1, 0.0));
    
    for (int i = 0; i < m; ++i){
        int idx = i + 1;

        double A = h[idx - 1] / 6.0;
        double B = (h[idx - 1] + h[idx]) / 3.0;
        double C = h[idx] / 6.0;
        double D = (y[idx + 1] - y[idx]) / h[idx] - (y[idx] - y[idx - 1]) / h[idx - 1];
        
        if (i > 0) system_matrix[i][i - 1] = A;
        system_matrix[i][i] = B;
        if (i < m - 1) system_matrix[i][i + 1] = C;
        system_matrix[i][m] = D;
    }

    std::vector<double> S = solve_tridiagonal_matrix(system_matrix, n);

    for (int i = 0; i < n-1; ++i){
        double a = y[i];
        double c = S[i] / 2.0;
        double d = (S[i + 1] - S[i]) / (6.0 * h[i]);
        double b = (y[i + 1] - y[i]) / h[i] - h[i]*(2.0 * S[i] + S[i+1]) / 6.0;
        f << a << ' ' << b << ' ' << c << ' ' << d << '\n';
    }
}

void lagrange_interpolation(std::vector<double> x, std::vector<double> y, std::vector<double> my_x, std::ofstream& f){
    int n = x.size();

    for (int k = 0; k < my_x.size(); ++k){
        double res = 0.0;
        for (int i = 0; i < n; ++i){
            double prod = y[i];
            for (int j = 0; j < n; ++j){
                if (j != i){
                    prod *= (my_x[k] - x[j]) / (x[i] - x[j]);
                }
            }
            res += prod;
        }

        f << res << '\n';
    }
}

int main(){
    int n;
    std::cout << "Введите количество узлов интерполирования: ";
    std::cin >> n;

    std::ofstream file_spline("spline.txt");
    std::ofstream file_func("func.txt");
    std::ofstream file_middles("mids.txt");
    std::ofstream file_lagrange("lagrange.txt");

    double a = 2.0;
    double b = 4.0;

    std::vector<double> nodes(n);
    std::vector<double> arr_y(n);
    std::vector<double> arr_h(n-1);
    std::vector<double> middles_h(n-1);
    std::vector<double> middles_h_y(n-1);

    for (int i = 0; i < n; ++i){
        double x = a + i *(b-a) / (n-1);
        nodes[i] = x;
        arr_y[i] = f(x);

        file_func << nodes[i] << ' ' << arr_y[i] << '\n';
    }

    for (int i = 0; i < n-1; ++i){
        arr_h[i] = nodes[i+1] - nodes[i];
        middles_h[i] = (nodes[i] + nodes[i+1]) / 2.0;
        middles_h_y[i] = f(middles_h[i]);

        file_middles << middles_h[i] << ' ' << middles_h_y[i] << '\n';
    }
         
    spline_interpolation(arr_h, arr_y, n, file_spline);     //n-1 уравнений, в каждом 4 неизвестных

    lagrange_interpolation(nodes, arr_y, middles_h, file_lagrange);      //Значения многочлена Лагранжа в серединах отрезков

    file_spline.close();
    file_lagrange.close();
    file_func.close();
    file_middles.close();

    return 0;
}