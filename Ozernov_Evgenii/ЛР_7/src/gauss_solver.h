#pragma once

#include <cstddef>
#include <string>
#include <vector>

// Стратегия выбора ведущего элемента.
// В этой лабораторной используется частичный выбор:
// в текущем столбце выбирается элемент с максимальным модулем.
enum class PivotStrategy {
    PartialColumnMax
};

// Тип результата решения СЛАУ.
enum class SolutionType {
    Unique,      // Единственное решение.
    Infinite,    // Бесконечно много решений.
    NoSolution,  // Решений нет (система несовместна).
    InvalidInput // Некорректный вход: размеры, структура матрицы и т.д.
};

// Итог работы метода Гаусса.
struct GaussResult {
    SolutionType type = SolutionType::InvalidInput; // Классификация результата.
    std::vector<double> x;                          // Найденный вектор решения (если применимо).
    std::size_t rank = 0;                           // Ранг матрицы коэффициентов.
    bool nonsingular = false;                       // Признак невырожденности (только для квадратной матрицы).
    std::string message;                            // Читаемое пояснение результата.
};

// Решение СЛАУ A * x = b методом Гаусса.
// Важно:
// 1) A и b передаются по значению, чтобы функция могла изменять их "на месте"
//    без создания отдельной расширенной матрицы [A|b].
//    Для больших задач желательно вызывать как solveGaussian(std::move(A), std::move(b), ...),
//    чтобы избежать лишнего копирования.
// 2) max_dimension задает рабочий предел по размерности задачи (по условию 10000).
GaussResult solveGaussian(
    std::vector<std::vector<double>> A,
    std::vector<double> b,
    PivotStrategy pivot = PivotStrategy::PartialColumnMax,
    double eps = 1e-12,
    std::size_t max_dimension = 10000
);
