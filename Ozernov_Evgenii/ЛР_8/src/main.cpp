#include "jacobi_solver.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

// Преобразование типа решения в читаемую строку.
std::string solutionTypeToString(SolutionType t) {
    switch (t) {
        case SolutionType::UNIQUE:
            return "Единственное решение";
        case SolutionType::INFINITE:
            return "Бесконечно много решений";
        case SolutionType::NONE:
            return "Решений нет";
        case SolutionType::UNDETERMINED:
            return "Не удалось строго классифицировать";
    }
    return "Неизвестный тип";
}

} // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Ниже простой интерфейс ввода:
    // 1) n
    // 2) матрица A (n x n)
    // 3) вектор b (n)
    // 4) стратегия выбора ведущего элемента: 0 или 1
    //
    // Пример:
    // 3
    // 10 1 1
    // 2 10 1
    // 2 2 10
    // 12 13 14
    // 1
    std::size_t n = 0;
    if (!(std::cin >> n)) {
        std::cerr << "Ошибка: не удалось прочитать размер матрицы n.\n";
        return 1;
    }

    SparseMatrix a(n);
    std::vector<double> b(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double value = 0.0;
            if (!(std::cin >> value)) {
                std::cerr << "Ошибка: не удалось прочитать элемент A[" << i << "][" << j << "].\n";
                return 1;
            }
            a.set(i, j, value);
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!(std::cin >> b[i])) {
            std::cerr << "Ошибка: не удалось прочитать элемент b[" << i << "].\n";
            return 1;
        }
    }

    int pivotChoice = 1;
    if (!(std::cin >> pivotChoice)) {
        std::cerr << "Ошибка: не удалось прочитать код стратегии ведущего элемента.\n";
        return 1;
    }

    JacobiConfig cfg;
    cfg.tolerance = 1e-10;
    cfg.maxIterations = 200000;
    cfg.zeroEps = 1e-12;
    cfg.exactRankCheckLimit = 250;
    cfg.pivotStrategy = (pivotChoice == 1)
                            ? PivotStrategy::PARTIAL_COLUMN_PIVOTING
                            : PivotStrategy::NONE;

    const SolveResult res = solveByJacobi(a, b, cfg);

    std::cout << "Статус: " << (res.success ? "Успех" : "Ошибка/Не сошёлся") << "\n";
    std::cout << "Сходимость: " << (res.converged ? "Да" : "Нет") << "\n";
    std::cout << "Тип решения: " << solutionTypeToString(res.solutionType) << "\n";
    std::cout << "Итераций: " << res.iterations << "\n";
    std::cout << "Достигнутая норма ||x(k+1)-x(k)||inf: " << res.achievedDelta << "\n";
    std::cout << "Оценка памяти (байт): " << res.estimatedMemoryBytes << "\n";
    std::cout << "Сообщение: " << res.message << "\n";

    if (res.success) {
        std::cout << "Решение x:\n";
        for (std::size_t i = 0; i < res.x.size(); ++i) {
            std::cout << "x[" << i << "] = " << res.x[i] << "\n";
        }
    }

    return 0;
}

