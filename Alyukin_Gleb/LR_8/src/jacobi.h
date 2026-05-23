#pragma once
#include <vector>

enum class SolveStatus {
    Unique,     // Сходимость достигнута
    NoSolution, // Расходится или деление на 0
    MaxIter     // Достигнуто макс. кол-во итераций без нужной точности
};

struct JacobiResult {
    SolveStatus status;
    std::vector<double> x;
    int iterations;
    double final_error;
};

// Структура для хранения элемента разреженной матрицы
struct SparseEntry {
    int col;
    double val;
};

JacobiResult jacobiSolve(std::vector<std::vector<double>> A, std::vector<double> b, 
                         int n, double eps = 1e-7, int max_iter = 1000);
