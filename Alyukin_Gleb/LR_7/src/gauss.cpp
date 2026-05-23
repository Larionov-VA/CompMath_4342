// g++ -O2 -std=c++17 gauss.cpp -o gauss
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>   // Для генерации случайных чисел
#include "gauss.h"

static constexpr double eps = GAUSS_EPS;

// Метод Гаусса с частичным выбором ведущего элемента
// A - матрица коэффициентов, размер m*n
// b - вектор правой части, размер m
// m - число уравнений
// n - число неизвестных
// Метод Гаусса с частичным выбором ведущего элемента
GaussResult gaussSolve(std::vector<std::vector<double>> A, std::vector<double> b, int m, int n)
{
    GaussResult res;
    res.rank = 0;
    res.nonsingular = false;

    std::vector<int> pivot_row(n, -1);
    int row = 0;

    // ПРЯМОЙ ХОД
    for (int col = 0; col < n && row < m; ++col) {
        int max_row = row;
        double max_val = std::abs(A[row][col]);

        for (int i = row + 1; i < m; ++i) {
            double val = std::abs(A[i][col]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val <= eps) continue;

        if (max_row != row) {
            std::swap(A[row], A[max_row]);
            std::swap(b[row], b[max_row]);
        }

        pivot_row[col] = row;

        for (int i = row + 1; i < m; ++i) {
            if (std::abs(A[i][col]) <= eps) {
                A[i][col] = 0.0;
                continue;
            }
            double factor = A[i][col] / A[row][col];
            for (int j = col + 1; j < n; ++j) {
                A[i][j] -= factor * A[row][j];
            }
            b[i] -= factor * b[row];
            A[i][col] = 0.0;
        }
        ++res.rank;
        ++row;
    }

    // ПРОВЕРКА НА НЕСОВМЕСТНОСТЬ
    for (int i = row; i < m; ++i) {
        bool row_is_zero = true;
        for (int j = 0; j < n; ++j) {
            if (std::abs(A[i][j]) > eps) { row_is_zero = false; break; }
        }
        if (row_is_zero && std::abs(b[i]) > eps) {
            res.status      = SolveStatus::NoSolution;
            res.nonsingular = false;
            return res;
        }
    }

    // ОБРАТНЫЙ ХОД
    res.x.assign(n, 0.0);
    for (int col = n - 1; col >= 0; --col) {
        if (pivot_row[col] == -1) continue;
        int cur = pivot_row[col];
        double sum = b[cur];
        for (int j = col + 1; j < n; ++j) {
            sum -= A[cur][j] * res.x[j];
        }
        res.x[col] = sum / A[cur][col];
    }

    // ИТОГОВАЯ КЛАССИФИКАЦИЯ
    if (res.rank == n && m == n) {
        res.status      = SolveStatus::Unique;
        res.nonsingular = true;
    } else if (res.rank == n) {
        res.status      = SolveStatus::Unique;
        res.nonsingular = false;
    } else {
        res.status      = SolveStatus::Infinite;
        res.nonsingular = false;
    }

    return res;
}

int main() {
    int m, n;
    std::cout << "Количество уравнений (m) и неизвестных (n): ";
    if (!(std::cin >> m >> n) || m <= 0 || n <= 0) {
        std::cout << "Ошибка: некорректные размеры матрицы\n";
        return 1;
    }

    std::vector<std::vector<double>> A(m, std::vector<double>(n));
    std::vector<double> b(m);

    std::cout << "Выберите способ заполнения:\n";
    std::cout << "1) Сгенерировать случайные элементы\n";
    std::cout << "2) Ввести вручную\n";
    std::cout << "Выбор: ";
    int choice;
    std::cin >> choice;

    if (choice == 1) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-10.0, 10.0);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = dist(gen);
            }
            b[i] = dist(gen);
        }

        std::cout << "\nСгенерированная матрица [A | b]:\n";
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << std::setw(10) << std::fixed << std::setprecision(3) << A[i][j];
            }
            std::cout << " | " << std::setw(10) << b[i] << "\n";
        }
        std::cout << std::endl;

    } else if (choice == 2) {
        std::cout << "Введите элементы матрицы A:\n";
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                if (!(std::cin >> A[i][j])) return 1;

        std::cout << "Введите элементы вектора b:\n";
        for (int i = 0; i < m; ++i)
            if (!(std::cin >> b[i])) return 1;
    } else {
        std::cout << "Неверный выбор.\n";
        return 1;
    }

    GaussResult res = gaussSolve(A, b, m, n);

    switch (res.status) {
        case SolveStatus::Unique:
            std::cout << "статус: единственное решение\n";  break;
        case SolveStatus::Infinite:
            std::cout << "статус: бесконечно много решений\n";  break;
        case SolveStatus::NoSolution:
            std::cout << "статус: нет решений\n";   break;
    }

    std::cout << "ранг: "          << res.rank                          << "\n";
    std::cout << "невырожденная: " << (res.nonsingular ? "да" : "нет")  << "\n";

    if (res.status == SolveStatus::Unique || res.status == SolveStatus::Infinite) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "вектор x:";
        for (double xi : res.x) std::cout << " " << xi;
        std::cout << "\n";

        if (res.status == SolveStatus::Infinite)
            std::cout << "(частное решение: свободные переменные приняты за 0)\n";
    }

    return 0;
}