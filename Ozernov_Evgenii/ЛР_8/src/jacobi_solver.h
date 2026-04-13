#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

// Перечисление стратегий выбора ведущего элемента.
// В текущей реализации реализованы:
// 1) NONE - перестановки строк не делаются.
// 2) PARTIAL_COLUMN_PIVOTING - для каждого столбца i выбирается строка с
//    максимальным по модулю элементом в этом столбце и переставляется наверх.
enum class PivotStrategy {
    NONE = 0,
    PARTIAL_COLUMN_PIVOTING = 1
};

// Разреженная матрица в формате "список строк".
// Для каждой строки храним пары (столбец, значение) только для ненулевых
// коэффициентов. Это позволяет комфортно работать с размерами 1000 и 10000.
class SparseMatrix {
public:
    explicit SparseMatrix(std::size_t n = 0);

    std::size_t size() const;

    // Установка значения A[row][col].
    // Если |value| очень мало, элемент удаляется из разреженного хранения.
    void set(std::size_t row, std::size_t col, double value);

    // Получение значения A[row][col]. Если элемента в хранилище нет, это 0.
    double get(std::size_t row, std::size_t col) const;

    // Доступ к списку ненулевых элементов конкретной строки.
    const std::vector<std::pair<std::size_t, double>>& row(std::size_t row) const;

    // Перестановка двух строк (нужна для выбора ведущего элемента).
    void swapRows(std::size_t r1, std::size_t r2);

    // Оценка памяти, занимаемой матрицей, в байтах.
    std::size_t memoryBytes() const;

private:
    std::vector<std::vector<std::pair<std::size_t, double>>> rows_;
};

// Классификация системы по наличию решений.
enum class SolutionType {
    UNIQUE,        // Единственное решение
    INFINITE,      // Бесконечно много решений
    NONE,          // Решений нет
    UNDETERMINED   // Однозначно определить по быстрой проверке не удалось
};

// Параметры метода Якоби.
struct JacobiConfig {
    double tolerance = 1e-10;                  // Критерий остановки по норме ||x(k+1)-x(k)||inf
    std::size_t maxIterations = 100000;        // Максимум итераций
    PivotStrategy pivotStrategy = PivotStrategy::PARTIAL_COLUMN_PIVOTING;
    double zeroEps = 1e-12;                    // Порог "считаем нулём"
    std::size_t exactRankCheckLimit = 250;     // До какого n делаем точную ранговую проверку
};

// Результат работы решателя.
struct SolveResult {
    bool success = false;                      // Удалось ли получить решение
    bool converged = false;                    // Сошёлся ли итерационный процесс
    SolutionType solutionType = SolutionType::UNIQUE;
    std::size_t iterations = 0;                // Количество выполненных итераций
    std::vector<double> x;                     // Найденный вектор решения (если есть)
    std::string message;                       // Подробный статус
    double achievedDelta = 0.0;                // Фактическая ||x(k+1)-x(k)||inf на последней итерации
    std::size_t estimatedMemoryBytes = 0;      // Оценка использованной памяти
};

// Основная функция решения СЛУ Ax=b методом Якоби.
SolveResult solveByJacobi(SparseMatrix a, std::vector<double> b, const JacobiConfig& config);

// Вспомогательная функция для генерации матрицы Гильберта размера n.
SparseMatrix buildHilbertMatrix(std::size_t n);

// Вспомогательная функция для генерации трёхдиагональной матрицы:
// diagValue на главной диагонали, offDiagValue на соседних.
SparseMatrix buildTridiagonalMatrix(std::size_t n, double diagValue, double offDiagValue);

