#pragma once
#include <vector>

inline constexpr double GAUSS_EPS = 1e-12;

enum class SolveStatus {
    Unique,       // единственное решение
    Infinite,     // бесконечно много решений (есть свободные переменные)
    NoSolution,   // система несовместна (противоречие 0 = c ≠ 0)
    InputError    // некорректные входные данные
};

struct GaussResult {
    SolveStatus         status;       // тип решения
    int                 rank;         // ранг матрицы A
    bool                nonsingular;  // true, если m==n и rank==n
    std::vector<double> x;            // вектор решения; при Infinite — частное
                                      // решение (свободные переменные = 0);
                                      // пуст при NoSolution / InputError
};

// A - матрица коэффициентов, размер m*n
// b - вектор правой части, размер m
// m - число уравнений
// n - число неизвестных
GaussResult gaussSolve(std::vector<std::vector<double>> A,
                       std::vector<double>              b,
                       int m, int n);