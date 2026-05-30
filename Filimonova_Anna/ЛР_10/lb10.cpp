#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

long double f(long double x) {
    //return std::exp(-x) - x*x*x;
    return 1/(1+25*x*x);
}

std::vector<long double> linspace(long double a, long double b, int N) {
    std::vector<long double> res(N);
    long double h = (b - a) / (N - 1);
    for (int i = 0; i < N; ++i)
        res[i] = a + i * h;
    return res;
}

std::vector<long double> divided_diffs(const std::vector<long double>& x,
                                       const std::vector<long double>& y) {
    int n = x.size();
    std::vector<long double> coeffs = y;
    for (int i = 1; i < n; ++i) {
        for (int j = n - 1; j >= i; --j) {
            coeffs[j] = (coeffs[j] - coeffs[j-1]) / (x[j] - x[j-i]);
        }
    }
    return coeffs;
}

long double newton_eval(long double t, const std::vector<long double>& x,
                        const std::vector<long double>& coeffs) {
    int n = coeffs.size();
    long double res = coeffs[0];
    long double term = 1.0L;
    for (int k = 1; k < n; ++k) {
        term *= (t - x[k-1]);
        res  += coeffs[k] * term;
    }
    return res;
}

std::vector<double> solve_tridiag(std::vector<double> lower,
                                  std::vector<double> diag,
                                  std::vector<double> upper,
                                  std::vector<double> rhs) {
    int n = diag.size();
    for (int i = 1; i < n; ++i) {
        double w = lower[i] / diag[i-1];
        diag[i] -= w * upper[i-1];
        rhs[i]  -= w * rhs[i-1];
    }
    std::vector<double> x(n);
    x[n-1] = rhs[n-1] / diag[n-1];
    for (int i = n-2; i >= 0; --i)
        x[i] = (rhs[i] - upper[i] * x[i+1]) / diag[i];
    return x;
}

struct SplineCoeffs {
    std::vector<double> a, b, c, d;
};

SplineCoeffs build_spline(const std::vector<double>& x,
                          const std::vector<double>& y) {
    int n = x.size() - 1;
    std::vector<double> h(n);
    for (int i = 0; i < n; ++i) h[i] = x[i+1] - x[i];

    std::vector<double> lower(n+1, 0.0), diag(n+1), upper(n+1, 0.0), rhs(n+1, 0.0);
    diag[0] = 1.0;
    for (int i = 1; i < n; ++i) {
        lower[i] = h[i-1];
        diag[i]  = 2.0 * (h[i-1] + h[i]);
        upper[i] = h[i];
        rhs[i]   = 3.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
    }
    diag[n] = 1.0;

    std::vector<double> c = solve_tridiag(lower, diag, upper, rhs);

    SplineCoeffs sp;
    sp.a.resize(n); sp.b.resize(n); sp.c.resize(n); sp.d.resize(n);
    for (int i = 0; i < n; ++i) {
        sp.a[i] = y[i];
        sp.b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (2.0*c[i] + c[i+1]) / 3.0;
        sp.c[i] = c[i];
        sp.d[i] = (c[i+1] - c[i]) / (3.0 * h[i]);
    }
    return sp;
}

double spline_eval(double t, const std::vector<double>& x,
                   const SplineCoeffs& sp) {
    int n = x.size() - 1;
    for (int i = 0; i < n; ++i) {
        if (t >= x[i] && t <= x[i+1]) {
            double dx = t - x[i];
            return sp.a[i] + sp.b[i]*dx + sp.c[i]*dx*dx + sp.d[i]*dx*dx*dx;
        }
    }
    return 0.0;
}

int main() {
    long double a, b;
    std::cout << "Границы отрезка a и b: ";
    std::cin >> a >> b;
    if (a >= b) { std::cout << "a должно быть меньше b\n"; return 1; }

    std::cout << "Количество отрезков разбиения: ";
    int n;
    std::cin >> n;
    if (n < 2) { std::cout << "Число отрезков должно быть >= 2.\n"; return 1; }

    std::vector<long double> x_nodes = linspace(a, b, n+1);
    std::vector<long double> y_nodes(n+1);
    for (int i = 0; i <= n; ++i)
        y_nodes[i] = f(x_nodes[i]);

    std::vector<long double> coeffs = divided_diffs(x_nodes, y_nodes);

    std::vector<double> x_nodes_d(n+1), y_nodes_d(n+1);
    for (int i = 0; i <= n; ++i) {
        x_nodes_d[i] = (double)x_nodes[i];
        y_nodes_d[i] = (double)y_nodes[i];
    }
    SplineCoeffs sp = build_spline(x_nodes_d, y_nodes_d);

    int N_plot = 1001;
    std::vector<long double> x_plot = linspace(a, b, N_plot);
    long double max_err_poly = 0.0L, max_err_spline = 0.0L;

    std::ofstream plot("plot_data.dat");
    for (int i = 0; i < N_plot; ++i) {
        long double xi = x_plot[i];
        long double fi = f(xi);
        long double poly = newton_eval(xi, x_nodes, coeffs);
        double spl = spline_eval((double)xi, x_nodes_d, sp);

        plot << (double)xi << " " << (double)fi << " "
             << (double)poly << " " << spl << "\n";

        long double e1 = std::fabs(fi - poly);
        long double e2 = std::fabs(fi - (long double)spl);
        if (e1 > max_err_poly) max_err_poly = e1;
        if (e2 > max_err_spline) max_err_spline = e2;
    }
    plot.close();

    std::cout << std::fixed << std::setprecision(18);
    std::cout << "\nКоэффициенты для n = " << n << "\n\n";
    std::cout << "Разделённые разности полинома Ньютона:\n";
    for (int i = 0; i <= n; ++i)
        std::cout << "f[x0..x" << i << "] = " << coeffs[i] << "\n";

    std::cout << "\nКоэффициенты кубического сплайна (a,b,c,d) по сегментам:\n";
    for (int i = 0; i < n; ++i) {
        std::cout << "Сегмент " << i << ": ["
                  << x_nodes_d[i] << ", " << x_nodes_d[i+1] << "]\n"
                  << "  a = " << sp.a[i] << "  b = " << sp.b[i]
                  << "  c = " << sp.c[i] << "  d = " << sp.d[i] << "\n";
    }

    std::cout << "\n\nМаксимальная погрешность полинома: " << max_err_poly << "\n";
    std::cout << "Максимальная погрешность сплайна: " << max_err_spline << "\n";
    return 0;
}