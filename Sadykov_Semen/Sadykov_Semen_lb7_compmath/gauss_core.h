#ifndef GAUSS_CORE_H
#define GAUSS_CORE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <string>

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

struct MatrixStats {
    int rank;
    bool is_degenerate;
    double memory_bytes;
};

class GaussSolver {
public:
    // 1. Анализ матрицы (вычисление ранга методом Гаусса)
    MatrixStats analyzeMatrix(Matrix A) {
        int n = A.size();
        int m = A[0].size();
        int rank = 0;
        vector<bool> row_used(n, false);

        for (int i = 0; i < m; ++i) {
            int j;
            for (j = 0; j < n; ++j) {
                if (!row_used[j] && abs(A[j][i]) > 1e-12) break;
            }

            if (j != n) {
                rank++;
                row_used[j] = true;
                for (int k = 0; k < n; ++k) {
                    if (k != j && abs(A[k][i]) > 1e-12) {
                        double val = A[k][i] / A[j][i];
                        for (int p = i + 1; p < m; ++p)
                            A[k][p] -= A[j][p] * val;
                    }
                }
            }
        }
        double mem = (double)(n * m * sizeof(double) + n * sizeof(double));
        return {rank, (n == m ? rank < n : false), mem};
    }

    // 2. Классический метод Гаусса с выбором ведущего элемента
    Vector solve(Matrix A, Vector b, bool& success) {
        int n = A.size();
        int m = A[0].size();
        Vector x(m, 0.0);

        for (int i = 0; i < m; i++) {
            // Поиск ведущего элемента
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (abs(A[j][i]) > abs(A[pivot][i])) pivot = j; 
            }

            if (abs(A[pivot][i]) < 1e-18) {
                success = false;
                return x;
            }

            swap(A[i], A[pivot]);
            swap(b[i], b[pivot]);

            // Прямой ход
            for (int j = i + 1; j < n; j++) {
                double factor = A[j][i] / A[i][i];
                b[j] -= factor * b[i];
                for (int k = i; k < m; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }

        // Обратный ход
        for (int i = m - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < m; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        success = true;
        return x;
    }

    // 3. Метод Прогонки (для специализированных тестов)
    Vector solveTridiagonal(const vector<double>& a, const vector<double>& b, const vector<double>& c, vector<double> d, bool& success) {
        int n = d.size();
        Vector x(n);
        vector<double> c_prime(n), d_prime(n);

        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];

        for (int i = 1; i < n; i++) {
            double m = 1.0 / (b[i] - a[i-1] * c_prime[i-1]);
            if (i < n - 1) c_prime[i] = c[i] * m;
            d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) * m;
        }

        x[n-1] = d_prime[n-1];
        for (int i = n - 2; i >= 0; i--) {
            x[i] = d_prime[i] - c_prime[i] * x[i+1];
        }
        success = true;
        return x;
    }

    double getResidual(const Matrix& A, const Vector& x, const Vector& b) {
        double max_err = 0;
        for (int i = 0; i < A.size(); ++i) {
            double s = 0;
            for (int j = 0; j < x.size(); ++j) s += A[i][j] * x[j];
            max_err = max(max_err, abs(s - b[i]));
        }
        return max_err;
    }
};

#endif