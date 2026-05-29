#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

 
// Перечисление для типа решения системы.
 
enum class SolutionType {
    Unique,   // Единственное решение
    Infinite, // Бесконечно много решений
    None      // Решений нет
};

 
// Структура с результатом работы метода Гаусса.
 
struct GaussResult {
    SolutionType type = SolutionType::None;

    // Вектор решения.
    // Для Unique здесь будет единственное решение.
    // Для Infinite здесь можно хранить одно частное решение
    // (свободные переменные положены равными нулю).
    // Для None вектор будет пустым.
    std::vector<double> x;

    // Ранг матрицы коэффициентов A.
    int rankA = 0;

    // Ранг расширенной матрицы [A | b].
    int rankAugmented = 0;

    // Признак вырожденности матрицы A.
    // Если rankA < n, то матрица вырождена.
    bool isDegenerate = false;
};

 
// Класс, реализующий решение СЛУ методом Гаусса.
 
class GaussianSolver {
public:
    // Конструктор.
    // eps - допустимая погрешность для сравнения с нулём.
    explicit GaussianSolver(double eps = 1e-12) : eps_(eps) {}

     
    // Основной метод решения системы A * x = b.
    //
    // A - квадратная матрица размера n x n
    // b - вектор правых частей размера n
    //
    // Метод:
    // 1) Проверяет корректность входных данных
    // 2) Строит расширенную матрицу [A | b]
    // 3) Выполняет прямой ход Гаусса
    // 4) Во время прямого хода фактически проверяет
    //    невырожденность матрицы:
    //    если в очередном столбце не найден ведущий элемент,
    //    значит ранг матрицы меньше n
    // 5) После прямого хода определяет:
    //    - есть ли решение,
    //    - единственно ли оно,
    //    - или решений бесконечно много
    // 6) При необходимости делает обратный ход
     
    GaussResult solve(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b) const
    {
        validateInput(A, b);

        const int n = static_cast<int>(A.size());

        // Строим расширенную матрицу [A | b].
        std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                aug[i][j] = A[i][j];
            }
            aug[i][n] = b[i];
        }

        // row - номер текущей строки, в которую будем ставить
        // очередной ведущий элемент.
        int row = 0;

        // Прямой ход метода Гаусса.
        // Перебираем столбцы слева направо.
        for (int col = 0; col < n && row < n; ++col) {
            // Ищем строку с максимальным по модулю элементом
            // в текущем столбце среди строк row..n-1.
            int pivotRow = selectPivotRow(aug, row, col);

            // Если в этом столбце нет ненулевого элемента,
            // то ведущий элемент поставить нельзя.
            // Это означает, что столбец не даёт новый pivot,
            // а матрица потенциально вырождена.
            if (pivotRow == -1) {
                continue;
            }

            // Если нужно, меняем строки местами,
            // чтобы pivot оказался в строке row.
            if (pivotRow != row) {
                std::swap(aug[pivotRow], aug[row]);
            }

            // Обнуляем элементы ниже ведущего.
            eliminateBelow(aug, row, col);

            // Переходим к следующей строке.
            ++row;
        }

        // После прямого хода вычисляем ранги матрицы A
        // и расширенной матрицы [A | b].
        const auto [rankA, rankAugmented] = computeRanks(aug);

        GaussResult result;
        result.rankA = rankA;
        result.rankAugmented = rankAugmented;
        result.isDegenerate = (rankA < n);

        // Если rank(A) < rank([A|b]), то система несовместна.
        if (rankA < rankAugmented) {
            result.type = SolutionType::None;
            return result;
        }

        // Если rank(A) = rank([A|b]) < n, то решений бесконечно много.
        // В этом случае построим одно частное решение:
        // свободные переменные положим равными нулю.
        if (rankA < n) {
            result.type = SolutionType::Infinite;
            result.x = buildOneSolutionFromEchelon(aug);
            return result;
        }

        // Если rank(A) = rank([A|b]) = n, то решение единственно.
        result.type = SolutionType::Unique;
        result.x = buildOneSolutionFromEchelon(aug);
        return result;
    }

private:
    double eps_;

     
    // Проверка "почти ноль".
     
    bool isZero(double value) const {
        return std::fabs(value) < eps_;
    }

     
    // Проверка корректности входных данных.
    //
    // Требования:
    // - A должна быть квадратной
    // - размер b должен совпадать с числом строк A
     
    void validateInput(const std::vector<std::vector<double>>& A,
                       const std::vector<double>& b) const
    {
        if (A.empty()) {
            throw std::invalid_argument("Матрица A не должна быть пустой.");
        }

        const std::size_t n = A.size();

        if (b.size() != n) {
            throw std::invalid_argument("Размер вектора b не совпадает с размером матрицы A.");
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (A[i].size() != n) {
                throw std::invalid_argument("Матрица A должна быть квадратной.");
            }
        }
    }

     
    // Выбор ведущего элемента.
    //
    // Метод выбора:
    // частичный выбор по столбцу.
    //
    // Среди строк startRow..n-1 ищем строку с максимальным
    // по модулю элементом в столбце col.
    //
    // Если все элементы в этом столбце близки к нулю,
    // возвращаем -1.
     
    int selectPivotRow(const std::vector<std::vector<double>>& aug,
                       int startRow,
                       int col) const
    {
        const int n = static_cast<int>(aug.size());

        int bestRow = -1;
        double bestAbsValue = 0.0;

        for (int i = startRow; i < n; ++i) {
            const double currentAbsValue = std::fabs(aug[i][col]);
            if (currentAbsValue > bestAbsValue) {
                bestAbsValue = currentAbsValue;
                bestRow = i;
            }
        }

        if (bestRow == -1 || isZero(bestAbsValue)) {
            return -1;
        }

        return bestRow;
    }

     
    // Прямой ход: обнуление элементов ниже ведущего.
    //
    // Пусть pivot находится в позиции (pivotRow, pivotCol).
    // Тогда для каждой строки i > pivotRow:
    //
    // row_i = row_i - factor * row_pivot,
    // где factor = aug[i][pivotCol] / aug[pivotRow][pivotCol]
    //
    // Мы специально пропускаем строки, где элемент уже близок к нулю,
    // чтобы не делать лишнюю работу.
     
    void eliminateBelow(std::vector<std::vector<double>>& aug,
                        int pivotRow,
                        int pivotCol) const
    {
        const int n = static_cast<int>(aug.size());
        const double pivot = aug[pivotRow][pivotCol];

        for (int i = pivotRow + 1; i < n; ++i) {
            if (isZero(aug[i][pivotCol])) {
                // Если под pivot уже стоит 0,
                // то ничего делать не нужно.
                continue;
            }

            const double factor = aug[i][pivotCol] / pivot;

            // После вычитания этот элемент должен стать нулём.
            aug[i][pivotCol] = 0.0;

            // Изменяем элементы справа от pivot.
            for (int j = pivotCol + 1; j <= n; ++j) {
                aug[i][j] -= factor * aug[pivotRow][j];

                // Малые вычислительные шумы прижимаем к нулю,
                // чтобы матрица оставалась "чище".
                if (isZero(aug[i][j])) {
                    aug[i][j] = 0.0;
                }
            }
        }
    }

     
    // Вычисление рангов матрицы A и расширенной матрицы [A | b]
    // по уже приведённой ступенчатой форме.
    //
    // Для каждой строки:
    // - если есть ненулевой коэффициент среди столбцов 0..n-1,
    //   то эта строка увеличивает rank(A)
    // - если в строке есть ненулевой коэффициент среди 0..n,
    //   то она увеличивает rank([A|b])
    //
    // Особенно важен случай:
    // [0 0 0 ... 0 | c], где c != 0
    // Тогда rank([A|b]) > rank(A), и система несовместна.
     
    std::pair<int, int> computeRanks(const std::vector<std::vector<double>>& aug) const
    {
        const int n = static_cast<int>(aug.size());

        int rankA = 0;
        int rankAugmented = 0;

        for (int i = 0; i < n; ++i) {
            bool hasNonZeroInA = false;
            bool hasNonZeroInAugmented = false;

            for (int j = 0; j < n; ++j) {
                if (!isZero(aug[i][j])) {
                    hasNonZeroInA = true;
                    hasNonZeroInAugmented = true;
                    break;
                }
            }

            // Если в коэффициентах всё нули, проверяем свободный член.
            if (!hasNonZeroInA && !isZero(aug[i][n])) {
                hasNonZeroInAugmented = true;
            }

            if (hasNonZeroInA) {
                ++rankA;
            }
            if (hasNonZeroInAugmented) {
                ++rankAugmented;
            }
        }

        return {rankA, rankAugmented};
    }

     
    // Находит индекс первого ненулевого коэффициента в строке
    // среди столбцов 0..n-1.
    //
    // Если таких нет, возвращает -1.
     
    int firstNonZeroColumn(const std::vector<double>& row) const
    {
        const int n = static_cast<int>(row.size()) - 1; // последний столбец - это b
        for (int j = 0; j < n; ++j) {
            if (!isZero(row[j])) {
                return j;
            }
        }
        return -1;
    }

     
    // Построение одного решения из ступенчатой формы.
    //
    // Этот метод работает и для случая единственного решения,
    // и для случая бесконечного числа решений.
    //
    // Идея:
    // - свободные переменные оставляем равными нулю
    // - для ведущих переменных идём снизу вверх
    //   и выполняем обратную подстановку
     
    std::vector<double> buildOneSolutionFromEchelon(
        const std::vector<std::vector<double>>& aug) const
    {
        const int n = static_cast<int>(aug.size());
        std::vector<double> x(n, 0.0);

        // Идём снизу вверх по строкам.
        for (int i = n - 1; i >= 0; --i) {
            const int leadCol = firstNonZeroColumn(aug[i]);

            // Полностью нулевая строка пропускается.
            if (leadCol == -1) {
                continue;
            }

            // Считаем сумму уже найденных членов справа.
            double sum = aug[i][n];
            for (int j = leadCol + 1; j < n; ++j) {
                sum -= aug[i][j] * x[j];
            }

            x[leadCol] = sum / aug[i][leadCol];
        }

        return x;
    }
};

 
// Вспомогательная функция печати решения.
// Не обязательна для алгоритма, но удобна в main.cpp.
 
inline void printGaussResult(const GaussResult& result)
{
    std::cout << std::setprecision(15);

    std::cout << "rank(A) = " << result.rankA << "\n";
    std::cout << "rank([A|b]) = " << result.rankAugmented << "\n";
    std::cout << "Матрица A "
              << (result.isDegenerate ? "вырождена" : "невырождена")
              << "\n";

    if (result.type == SolutionType::None) {
        std::cout << "Система несовместна: решений нет.\n";
        return;
    }

    if (result.type == SolutionType::Infinite) {
        std::cout << "Система имеет бесконечно много решений.\n";
        std::cout << "Одно частное решение (свободные переменные = 0):\n";
    } else {
        std::cout << "Система имеет единственное решение:\n";
    }

    for (std::size_t i = 0; i < result.x.size(); ++i) {
        std::cout << "x[" << i << "] = " << result.x[i] << "\n";
    }
}