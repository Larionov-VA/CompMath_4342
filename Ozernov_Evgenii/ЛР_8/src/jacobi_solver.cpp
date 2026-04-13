#include "jacobi_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

// Служебная структура с результатом рангового анализа.
struct RankAnalysisResult {
    SolutionType type = SolutionType::UNDETERMINED;
    std::string message;
};

// Проверка: есть ли в строке только нули.
bool isZeroRow(const std::vector<std::pair<std::size_t, double>>& row, double eps) {
    for (const auto& item : row) {
        if (std::abs(item.second) > eps) {
            return false;
        }
    }
    return true;
}

// Быстрый контроль на противоречивость системы:
// если есть строка 0 ... 0 | bi, где bi != 0, решений точно нет.
bool hasInconsistentZeroRow(const SparseMatrix& a, const std::vector<double>& b, double eps) {
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (isZeroRow(a.row(i), eps) && std::abs(b[i]) > eps) {
            return true;
        }
    }
    return false;
}

// Проверка строгого диагонального преобладания по строкам:
// |a_ii| > sum_{j!=i} |a_ij| для каждой строки i.
// Для таких матриц:
// 1) матрица невырождена,
// 2) метод Якоби гарантированно сходится.
bool isStrictlyDiagonallyDominant(const SparseMatrix& a, double eps) {
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        double diagAbs = 0.0;
        double offDiagSum = 0.0;
        for (const auto& [col, value] : a.row(i)) {
            if (col == i) {
                diagAbs = std::abs(value);
            } else {
                offDiagSum += std::abs(value);
            }
        }
        if (!(diagAbs > offDiagSum + eps)) {
            return false;
        }
    }
    return true;
}

// Перестановка строк по стратегии частичного выбора главного элемента по столбцу.
// Для каждого столбца i выбирается строка r >= i с максимальным |a_ri|.
void applyPartialColumnPivoting(SparseMatrix& a, std::vector<double>& b, double eps) {
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t pivotRow = i;
        double pivotAbs = std::abs(a.get(i, i));

        for (std::size_t r = i + 1; r < n; ++r) {
            const double candidate = std::abs(a.get(r, i));
            if (candidate > pivotAbs) {
                pivotAbs = candidate;
                pivotRow = r;
            }
        }

        // Если в столбце подходящий элемент не найден, строку не двигаем.
        // Ситуация потенциальной вырожденности будет дополнительно проверена позже.
        if (pivotAbs <= eps) {
            continue;
        }

        if (pivotRow != i) {
            a.swapRows(i, pivotRow);
            std::swap(b[i], b[pivotRow]);
        }
    }
}

// Точная проверка типа решений через ранги rank(A) и rank([A|b]).
// Используется только для относительно небольших размеров, потому что алгоритм O(n^3).
RankAnalysisResult analyzeByRanks(const SparseMatrix& a, const std::vector<double>& b, double eps) {
    const std::size_t n = a.size();
    std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1, 0.0));

    // Переносим разреженную матрицу в плотный вид для удобства гауссова исключения.
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& [col, value] : a.row(i)) {
            aug[i][col] = value;
        }
        aug[i][n] = b[i];
    }

    std::size_t row = 0;
    for (std::size_t col = 0; col < n && row < n; ++col) {
        std::size_t pivot = row;
        double maxAbs = std::abs(aug[row][col]);
        for (std::size_t r = row + 1; r < n; ++r) {
            const double curAbs = std::abs(aug[r][col]);
            if (curAbs > maxAbs) {
                maxAbs = curAbs;
                pivot = r;
            }
        }

        if (maxAbs <= eps) {
            continue;
        }

        if (pivot != row) {
            std::swap(aug[pivot], aug[row]);
        }

        for (std::size_t r = row + 1; r < n; ++r) {
            if (std::abs(aug[r][col]) <= eps) {
                continue;
            }
            const double factor = aug[r][col] / aug[row][col];
            for (std::size_t c = col; c <= n; ++c) {
                aug[r][c] -= factor * aug[row][c];
            }
        }

        ++row;
    }

    std::size_t rankA = 0;
    std::size_t rankAug = 0;
    for (std::size_t i = 0; i < n; ++i) {
        bool nonZeroA = false;
        for (std::size_t j = 0; j < n; ++j) {
            if (std::abs(aug[i][j]) > eps) {
                nonZeroA = true;
                break;
            }
        }

        const bool nonZeroB = std::abs(aug[i][n]) > eps;
        if (nonZeroA) {
            ++rankA;
        }
        if (nonZeroA || nonZeroB) {
            ++rankAug;
        }
    }

    RankAnalysisResult result;
    if (rankA < rankAug) {
        result.type = SolutionType::NONE;
        result.message = "Система несовместна: rank(A) < rank([A|b]).";
    } else if (rankA == rankAug && rankA < n) {
        result.type = SolutionType::INFINITE;
        result.message = "Система совместна, но имеет бесконечно много решений: rank(A)=rank([A|b])<n.";
    } else if (rankA == n) {
        result.type = SolutionType::UNIQUE;
        result.message = "Система имеет единственное решение: rank(A)=rank([A|b])=n.";
    } else {
        result.type = SolutionType::UNDETERMINED;
        result.message = "Не удалось классифицировать систему по рангу.";
    }
    return result;
}

// Форматирование вектора с ограничением количества выводимых элементов.
std::string vectorToString(const std::vector<double>& v, std::size_t maxElements = 8) {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < v.size() && i < maxElements; ++i) {
        if (i > 0) {
            oss << ", ";
        }
        oss << v[i];
    }
    if (v.size() > maxElements) {
        oss << ", ...";
    }
    oss << "]";
    return oss.str();
}

} // namespace

SparseMatrix::SparseMatrix(std::size_t n) : rows_(n) {}

std::size_t SparseMatrix::size() const {
    return rows_.size();
}

void SparseMatrix::set(std::size_t row, std::size_t col, double value) {
    if (row >= rows_.size() || col >= rows_.size()) {
        throw std::out_of_range("Индекс строки/столбца выходит за пределы матрицы.");
    }

    auto& r = rows_[row];
    auto it = std::find_if(r.begin(), r.end(), [col](const auto& item) { return item.first == col; });

    constexpr double storageEps = 1e-15;
    if (std::abs(value) <= storageEps) {
        if (it != r.end()) {
            r.erase(it);
        }
        return;
    }

    if (it != r.end()) {
        it->second = value;
    } else {
        r.emplace_back(col, value);
    }
}

double SparseMatrix::get(std::size_t row, std::size_t col) const {
    if (row >= rows_.size() || col >= rows_.size()) {
        throw std::out_of_range("Индекс строки/столбца выходит за пределы матрицы.");
    }
    const auto& r = rows_[row];
    auto it = std::find_if(r.begin(), r.end(), [col](const auto& item) { return item.first == col; });
    return (it == r.end()) ? 0.0 : it->second;
}

const std::vector<std::pair<std::size_t, double>>& SparseMatrix::row(std::size_t row) const {
    if (row >= rows_.size()) {
        throw std::out_of_range("Индекс строки выходит за пределы матрицы.");
    }
    return rows_[row];
}

void SparseMatrix::swapRows(std::size_t r1, std::size_t r2) {
    if (r1 >= rows_.size() || r2 >= rows_.size()) {
        throw std::out_of_range("Индекс строки выходит за пределы матрицы.");
    }
    if (r1 == r2) {
        return;
    }
    std::swap(rows_[r1], rows_[r2]);
}

std::size_t SparseMatrix::memoryBytes() const {
    std::size_t bytes = sizeof(*this);
    bytes += rows_.capacity() * sizeof(std::vector<std::pair<std::size_t, double>>);
    for (const auto& r : rows_) {
        bytes += r.capacity() * sizeof(std::pair<std::size_t, double>);
    }
    return bytes;
}

SolveResult solveByJacobi(SparseMatrix a, std::vector<double> b, const JacobiConfig& config) {
    SolveResult result;
    const std::size_t n = a.size();

    if (n == 0) {
        result.message = "Пустая матрица: решать нечего.";
        return result;
    }

    if (b.size() != n) {
        result.message = "Размер вектора b не совпадает с размером матрицы A.";
        return result;
    }

    // Применяем выбранную стратегию выбора ведущего элемента.
    if (config.pivotStrategy == PivotStrategy::PARTIAL_COLUMN_PIVOTING) {
        applyPartialColumnPivoting(a, b, config.zeroEps);
    }

    // Быстрый контроль несовместности.
    if (hasInconsistentZeroRow(a, b, config.zeroEps)) {
        result.solutionType = SolutionType::NONE;
        result.message = "Найдена строка вида 0 ... 0 | b_i (b_i != 0): решений нет.";
        return result;
    }

    const bool strictlyDominant = isStrictlyDiagonallyDominant(a, config.zeroEps);

    // Если строгая диагональная доминантность есть, решение единственно гарантированно.
    // Иначе для небольших n делаем точный ранговый анализ.
    if (strictlyDominant) {
        result.solutionType = SolutionType::UNIQUE;
    } else if (n <= config.exactRankCheckLimit) {
        const RankAnalysisResult rankResult = analyzeByRanks(a, b, config.zeroEps);
        result.solutionType = rankResult.type;
        if (rankResult.type == SolutionType::NONE || rankResult.type == SolutionType::INFINITE) {
            result.message = rankResult.message;
            return result;
        }
    } else {
        result.solutionType = SolutionType::UNDETERMINED;
    }

    // Начальное приближение x(0) = 0.
    std::vector<double> prevX(n, 0.0);
    std::vector<double> nextX(n, 0.0);

    // Оценка памяти на основные структуры данных, участвующие в расчёте.
    result.estimatedMemoryBytes =
        a.memoryBytes() +
        b.capacity() * sizeof(double) +
        prevX.capacity() * sizeof(double) +
        nextX.capacity() * sizeof(double);

    double previousDelta = std::numeric_limits<double>::infinity();

    for (std::size_t iter = 1; iter <= config.maxIterations; ++iter) {
        double maxDelta = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            double diag = 0.0;
            double sum = 0.0;

            // В формуле Якоби используем только значения из предыдущей итерации prevX.
            // x_i(k+1) = (b_i - sum_{j!=i}(a_ij*x_j(k))) / a_ii
            for (const auto& [col, value] : a.row(i)) {
                if (col == i) {
                    diag = value;
                } else {
                    sum += value * prevX[col];
                }
            }

            // Проверка невырожденности "в процессе работы":
            // если на диагонали ноль, формула Якоби для данной строки неприменима.
            if (std::abs(diag) <= config.zeroEps) {
                result.message =
                    "На итерации обнаружен нулевой (или близкий к нулю) диагональный элемент. "
                    "Матрица вырождена или некорректно переставлена для метода Якоби.";
                result.iterations = iter - 1;
                result.achievedDelta = maxDelta;
                return result;
            }

            const double newValue = (b[i] - sum) / diag;
            if (!std::isfinite(newValue)) {
                result.message = "Численная ошибка: получено нечисловое значение в итерациях.";
                result.iterations = iter - 1;
                result.achievedDelta = maxDelta;
                return result;
            }

            nextX[i] = newValue;
            maxDelta = std::max(maxDelta, std::abs(nextX[i] - prevX[i]));
        }

        result.iterations = iter;
        result.achievedDelta = maxDelta;

        if (maxDelta < config.tolerance) {
            result.success = true;
            result.converged = true;
            result.x = nextX;

            if (result.solutionType == SolutionType::UNDETERMINED) {
                result.message =
                    "Метод сошёлся, но для большого недиагонально-доминантного случая "
                    "тип решения не был строго классифицирован ранговой проверкой.";
            } else {
                result.message = "Метод Якоби успешно сошёлся. x = " + vectorToString(result.x);
            }
            return result;
        }

        // Простейший контроль на явное расхождение: если норма стала огромной, прерываем.
        if (iter > 10 && maxDelta > previousDelta * 1e6) {
            result.message =
                "Итерации расходятся (резкий рост нормы). Вероятно, для данной матрицы "
                "метод Якоби не подходит без дополнительного преобразования.";
            result.x = nextX;
            return result;
        }

        previousDelta = maxDelta;
        std::swap(prevX, nextX);
    }

    // До заданной точности не дошли, но частичный результат сохраним.
    result.success = false;
    result.converged = false;
    result.x = prevX;
    result.message =
        "Достигнут предел итераций без сходимости. Увеличьте maxIterations, ослабьте tolerance "
        "или преобразуйте систему (например, перестановками/масштабированием).";
    return result;
}

SparseMatrix buildHilbertMatrix(std::size_t n) {
    SparseMatrix h(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            h.set(i, j, 1.0 / static_cast<double>(i + j + 1));
        }
    }
    return h;
}

SparseMatrix buildTridiagonalMatrix(std::size_t n, double diagValue, double offDiagValue) {
    SparseMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) {
        m.set(i, i, diagValue);
        if (i > 0) {
            m.set(i, i - 1, offDiagValue);
        }
        if (i + 1 < n) {
            m.set(i, i + 1, offDiagValue);
        }
    }
    return m;
}

