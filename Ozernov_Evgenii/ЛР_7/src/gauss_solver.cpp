#include "gauss_solver.h"

#include <algorithm>
#include <cmath>

namespace {

// Проверка, что строка матрицы имеет вид [0 0 ... 0],
// то есть все коэффициенты по модулю меньше eps.
bool isZeroRow(const std::vector<double>& row, double eps) {
    for (double value : row) {
        if (std::fabs(value) > eps) {
            return false;
        }
    }
    return true;
}

} // namespace

GaussResult solveGaussian(
    std::vector<std::vector<double>> A,
    std::vector<double> b,
    PivotStrategy pivot,
    double eps,
    std::size_t max_dimension
) {
    GaussResult result;

    // 1. Валидация входа.

    // Матрица не должна быть пустой.
    if (A.empty() || A[0].empty()) {
        result.type = SolutionType::InvalidInput;
        result.message = "Матрица коэффициентов не должна быть пустой.";
        return result;
    }

    const std::size_t m = A.size();
    const std::size_t n = A[0].size();

    // Размер вектора правой части должен совпадать с числом строк матрицы коэффициентов.
    if (b.size() != m) {
        result.type = SolutionType::InvalidInput;
        result.message = "Размер вектора правой части должен совпадать с числом строк матрицы коэффициентов.";
        return result;
    }

    // Матрица должна быть прямоугольной.
    for (const auto& row : A) {
        if (row.size() != n) {
            result.type = SolutionType::InvalidInput;
            result.message = "Матрица коэффициентов должна быть прямоугольной (одинаковая длина всех строк).";
            return result;
        }
    }

    // По условию лабораторной поддерживаем системы до 10000x10000.
    // Здесь ограничение накладывается отдельно на число строк и столбцов.
    if (m > max_dimension || n > max_dimension) {
        result.type = SolutionType::InvalidInput;
        result.message = "Размерность системы превышает допустимый предел.";
        return result;
    }

    // where[col] = индекс строки, в которой стоит ведущий элемент для col,
    // либо -1, если в этом столбце ведущий элемент не найден.
    std::vector<int> where(n, -1);

    std::size_t row = 0;  // Текущая рабочая строка (позиция для следующего ведущего элемента).
    std::size_t rank = 0; // Ранг считаем прямо в процессе прямого хода.


    // 2. Прямой ход метода Гаусса (с выбором ведущего элемента).
  
    for (std::size_t col = 0; col < n && row < m; ++col) {
        // ---- 2.1 Выбор ведущего элемента ----
        // По стратегии PartialColumnMax ищем в текущем столбце элемент
        // с максимальным модулем среди строк row..m-1.
        std::size_t sel = row;
        if (pivot == PivotStrategy::PartialColumnMax) {
            for (std::size_t i = row + 1; i < m; ++i) {
                if (std::fabs(A[i][col]) > std::fabs(A[sel][col])) {
                    sel = i;
                }
            }
        }

        // Если в столбце нет значимого (по eps) элемента,
        // то ведущий элемент здесь отсутствует: столбец линейно зависим.
        // Это фиксирует вырожденность в процессе исключения.
        if (std::fabs(A[sel][col]) <= eps) {
            continue;
        }

        // ---- 2.2 Перестановка строк ----
        // Переносим выбранную строку на позицию текущей рабочей строки.
        std::swap(A[sel], A[row]);
        std::swap(b[sel], b[row]);
        where[col] = static_cast<int>(row);

        // ---- 2.3 Обнуление элементов ниже ведущего ----
        const double pivot_value = A[row][col];
        for (std::size_t i = row + 1; i < m; ++i) {
            const double factor = A[i][col] / pivot_value;

            // Если коэффициент практически нулевой, строка уже имеет ноль в col.
            if (std::fabs(factor) <= eps) {
                A[i][col] = 0.0;
                continue;
            }

            // Явно обнуляем элемент под ведущим.
            A[i][col] = 0.0;

            // Обновляем оставшуюся часть строки.
            for (std::size_t j = col + 1; j < n; ++j) {
                A[i][j] -= factor * A[row][j];
            }

            // Обновляем правую часть.
            b[i] -= factor * b[row];
        }

        // После нахождения очередного ведущего элемента
        // увеличиваем ранг и переходим к следующей строке.
        ++row;
        ++rank;
    }


    // 3. Проверка совместности после прямого хода.

    // Если встречается строка вида 0 ... 0 | c, где c != 0,
    // система несовместна и решений нет.
    for (std::size_t i = 0; i < m; ++i) {
        if (isZeroRow(A[i], eps) && std::fabs(b[i]) > eps) {
            result.type = SolutionType::NoSolution;
            result.rank = rank;
            result.nonsingular = false;
            result.message = "Система несовместна: обнаружена нулевая строка с ненулевой правой частью.";
            return result;
        }
    }

    // 4. Обратный ход (восстановление решения).

    std::vector<double> x(n, 0.0);

    // Идем по столбцам справа налево.
    for (int col = static_cast<int>(n) - 1; col >= 0; --col) {
        // Если в столбце нет ведущего элемента, это свободная переменная.
        // Для получения одного конкретного решения оставляем x[col] = 0.
        if (where[static_cast<std::size_t>(col)] == -1) {
            continue;
        }

        const std::size_t r = static_cast<std::size_t>(where[static_cast<std::size_t>(col)]);
        double sum = b[r];

        // Вычитаем уже найденные компоненты.
        for (std::size_t j = static_cast<std::size_t>(col) + 1; j < n; ++j) {
            sum -= A[r][j] * x[j];
        }

        x[static_cast<std::size_t>(col)] = sum / A[r][static_cast<std::size_t>(col)];
    }


    // 5. Финальная классификация результата.

    // Невырожденность: только для квадратной матрицы и только при rank == n.
    // Это условие проверяется по ходу работы через найденные ведущие элементы.
    const bool nonsingular = (m == n && rank == n);
    result.rank = rank;
    result.nonsingular = nonsingular;
    result.x = std::move(x);

    // Если ранг меньше числа неизвестных, есть свободные переменные => бесконечно много решений.
    if (rank < n) {
        result.type = SolutionType::Infinite;
        result.message = "Система совместна и имеет бесконечно много решений.";
        return result;
    }

    // Иначе решение единственное.
    result.type = SolutionType::Unique;
    result.message = "Система имеет единственное решение.";
    return result;
}
