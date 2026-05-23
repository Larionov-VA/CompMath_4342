#ifndef JACOBI_CORE_H
#define JACOBI_CORE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

// Структура для хранения результатов анализа матрицы
struct MatrixStats {
    int rank;
    bool is_degenerate;
    double memory_bytes;
};

class JacobiSolver {
public:
    double eps;
    int max_iter;

    JacobiSolver(double e = 1e-7, int iter = 10000) : eps(e), max_iter(iter) {}

    // 1. Метод вычисления ранга и проверки на вырожденность (методом Гаусса)
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
                for (int p = i + 1; p < m; ++p) A[j][p] /= A[j][i];
                for (int k = 0; k < n; ++k) {
                    if (k != j && abs(A[k][i]) > 1e-12) {
                        for (int p = i + 1; p < m; ++p)
                            A[k][p] -= A[j][p] * A[k][i];
                    }
                }
            }
        }
        double mem = (double)(n * m * sizeof(double) + n * sizeof(double));
        return {rank, (n == m ? rank < n : true), mem};
    }

    // 2. Метод выбора ведущего элемента (частичный выбор по столбцу)
    // Переставляет строки, чтобы на диагонали были максимальные элементы для сходимости
    void applyPivoting(Matrix& A, Vector& b) {
        int n = A.size();
        for (int i = 0; i < n; ++i) {
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (abs(A[k][i]) > abs(A[max_row][i])) max_row = k;
            }
            swap(A[i], A[max_row]);
            swap(b[i], b[max_row]);
        }
    }

    // 3. Основной алгоритм Якоби
    Vector solve(Matrix A, Vector b, bool& success) {
        int rows = A.size();
        int cols = A[0].size();
        
        // Для неквадратных систем (Тест 9) берем размер по количеству неизвестных
        int n = cols; 
        Vector x(n, 0.0);
        Vector x_new(n, 0.0);

        applyPivoting(A, b);

        for (int iter = 0; iter < max_iter; ++iter) {
            for (int i = 0; i < n; ++i) {
                // Проверка деления на ноль
                if (abs(A[i][i]) < 1e-18) {
                    success = false;
                    return x;
                }

                double sum = 0;
                for (int j = 0; j < n; ++j) {
                    if (i != j) sum += A[i][j] * x[j];
                }
                x_new[i] = (b[i] - sum) / A[i][i];
            }

            // Проверка сходимости по норме разности
            double diff = 0;
            for (int i = 0; i < n; ++i) diff += pow(x_new[i] - x[i], 2);
            
            if (sqrt(diff) < eps) {
                success = true;
                return x_new;
            }
            x = x_new;

            // Проверка на расходимость
            if (diff > 1e20) { success = false; return x; }
        }

        success = false;
        return x;
    }

    // Расчет невязки ||Ax - b||
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