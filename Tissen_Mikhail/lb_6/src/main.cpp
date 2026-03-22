#include "methods.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double LEFT = 3.0;
constexpr double RIGHT = 4.0;
constexpr double X0 = 3.5;
constexpr int NMAX = 1000;

// Для f(x) = exp(1 / x^2) - ln(x) на [3, 4]
// m = min |f'(x)|, M = max |f'(x)|, lambda = 2 / (m + M)
constexpr double M = 0.41611252361050843;
constexpr double m = 0.28326545184118310;
constexpr double LAMBDA = 2.0 / (m + M);

struct SpeedRow {
    double eps;
    int iters;
    double x;
    double absError;
    double residual;
};

struct ConditioningRow {
    double delta;
    int iters;
    double x;
    double absError;
    double residual;
    bool ok;
};

double f(double x) {
    if (x <= 0.0) {
        throw std::domain_error("f(x) is undefined for x <= 0");
    }
    return std::exp(1.0 / (x * x)) - std::log(x);
}

double fp(double x) {
    if (x <= 0.0) {
        throw std::domain_error("f'(x) is undefined for x <= 0");
    }
    return -2.0 * std::exp(1.0 / (x * x)) / (x * x * x) - 1.0 / x;
}

double phi(double x) {
    return x + LAMBDA * f(x);
}

double dphi(double x) {
    return 1.0 + LAMBDA * fp(x);
}

double phiRounded(double x, double delta) {
    return Round(phi(x), delta);
}

void WriteSpeedCsv(const std::string& filename,
                   const std::vector<SpeedRow>& rows) {
    std::ofstream out(filename);
    out << "eps,iters,x,abs_error,residual\n";
    out << std::setprecision(17);
    for (const auto& row : rows) {
        out << row.eps << ','
            << row.iters << ','
            << row.x << ','
            << row.absError << ','
            << row.residual << '\n';
    }
}

void WriteConditioningCsv(const std::string& filename,
                          const std::vector<ConditioningRow>& rows) {
    std::ofstream out(filename);
    out << "delta,status,iters,x,abs_error,residual\n";
    out << std::setprecision(17);
    for (const auto& row : rows) {
        out << row.delta << ','
            << (row.ok ? "ok" : "nmax") << ','
            << row.iters << ','
            << row.x << ','
            << row.absError << ','
            << row.residual << '\n';
    }
}

}  // локальное пространство имён

int main() {
    try {
        std::cout << std::fixed << std::setprecision(15);

        const double q = EstimateQ(dphi, LEFT, RIGHT);
        std::cout << "Function: f(x) = exp(1/x^2) - ln(x)\n";
        std::cout << "Interval: [" << LEFT << "; " << RIGHT << "]\n";
        std::cout << "Initial approximation x0 = " << X0 << "\n";
        std::cout << "lambda = " << LAMBDA << "\n";
        std::cout << "q = max |phi'(x)| on [" << LEFT << "; " << RIGHT << "] = " << q << "\n\n";

        const auto reference = ITER(phi, X0, 1e-14, NMAX);
        const double xref = reference.x;
        std::cout << "Reference run (delta = 0, eps = 1e-14):\n";
        std::cout << "  xref  = " << xref << "\n";
        std::cout << "  iters = " << reference.iters << "\n";
        std::cout << "  |f(xref)| = " << std::fabs(f(xref)) << "\n\n";

        const std::vector<double> epsValues = {
            1e-1, 1e-2, 1e-3, 1e-4,
            1e-5, 1e-6, 1e-7, 1e-8
        };

        std::vector<SpeedRow> speedRows;
        std::cout << "Speed of convergence:\n";
        for (double eps : epsValues) {
            const auto result = ITER(phi, X0, eps, NMAX);
            const SpeedRow row{
                eps,
                result.iters,
                result.x,
                std::fabs(result.x - xref),
                std::fabs(f(result.x))
            };
            speedRows.push_back(row);

            std::cout << "  eps = " << std::scientific << eps
                      << ", iters = " << std::fixed << row.iters
                      << ", x = " << std::setprecision(15) << row.x
                      << ", |x - xref| = " << std::scientific << row.absError
                      << ", |f(x)| = " << row.residual << '\n';
        }
        std::cout << '\n';

        const std::vector<double> deltaValues = {
            1e-1, 1e-2, 1e-3, 1e-4,
            1e-5, 1e-6, 1e-7, 1e-8
        };

        std::vector<ConditioningRow> conditioningRows;
        std::cout << "Conditioning experiment (eps = 1e-8):\n";
        for (double delta : deltaValues) {
            const auto result = ITER(
                [delta](double x) {
                    return phiRounded(x, delta);
                },
                X0,
                1e-8,
                NMAX
            );

            const ConditioningRow row{
                delta,
                result.iters,
                result.x,
                std::fabs(result.x - xref),
                std::fabs(f(result.x)),
                result.ok
            };
            conditioningRows.push_back(row);

            std::cout << "  delta = " << std::scientific << delta
                      << ", status = " << (row.ok ? "ok" : "nmax")
                      << ", iters = " << std::fixed << row.iters
                      << ", x = " << std::setprecision(15) << row.x
                      << ", |x - xref| = " << std::scientific << row.absError
                      << ", |f(x)| = " << row.residual << '\n';
        }

        WriteSpeedCsv("speed_results.csv", speedRows);
        WriteConditioningCsv("conditioning_results.csv", conditioningRows);

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
