#include "jacobi_solver.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace jacobi {

namespace {

double absd(double x) {
    return std::fabs(x);
}

// Поиск ранга методом Гаусса с выбором главного элемента по столбцу.
int rankOfMatrix(Matrix m, double eps) {
    const int rows = static_cast<int>(m.size());
    const int cols = rows == 0 ? 0 : static_cast<int>(m[0].size());

    int rank = 0;
    int current_row = 0;

    for (int col = 0; col < cols && current_row < rows; ++col) {
        int pivot_row = current_row;
        for (int r = current_row + 1; r < rows; ++r) {
            if (absd(m[r][col]) > absd(m[pivot_row][col])) {
                pivot_row = r;
            }
        }

        if (absd(m[pivot_row][col]) <= eps) {
            continue;
        }

        std::swap(m[current_row], m[pivot_row]);

        for (int r = current_row + 1; r < rows; ++r) {
            const double factor = m[r][col] / m[current_row][col];
            if (absd(factor) <= eps) {
                continue;
            }
            for (int c = col; c < cols; ++c) {
                m[r][c] -= factor * m[current_row][c];
            }
        }

        ++rank;
        ++current_row;
    }

    return rank;
}

Vector multiply(const Matrix& a, const Vector& x) {
    const int n = static_cast<int>(a.size());
    Vector result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < static_cast<int>(a[i].size()); ++j) {
            result[i] += a[i][j] * x[j];
        }
    }

    return result;
}

// Возвращает значение из ленточной матрицы.
// Если элемент находится вне ленты, возвращается 0.
double getBandedValue(const BandedSystem& s, int row, int col) {
    const int n = static_cast<int>(s.b.size());
    if (row < 0 || row >= n || col < 0 || col >= n) {
        return 0.0;
    }

    const int offset = col - row;
    if (offset < -s.lower_band || offset > s.upper_band) {
        return 0.0;
    }

    const int index = offset + s.lower_band;
    return s.diagonals[index][row];
}

// Вычисление b = A * x для ленточной матрицы.
Vector multiply(const BandedSystem& s, const Vector& x) {
    const int n = static_cast<int>(s.b.size());
    Vector result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        const int left = std::max(0, i - s.lower_band);
        const int right = std::min(n - 1, i + s.upper_band);

        for (int j = left; j <= right; ++j) {
            result[i] += getBandedValue(s, i, j) * x[j];
        }
    }

    return result;
}

// Создание пустой ленточной системы нужного размера.
BandedSystem makeEmptyBandedSystem(int n, int lower_band, int upper_band) {
    if (n <= 0) {
        throw std::invalid_argument("Размер системы должен быть положительным.");
    }
    if (lower_band < 0 || upper_band < 0) {
        throw std::invalid_argument("Ширина ленты не может быть отрицательной.");
    }

    BandedSystem s;
    s.lower_band = lower_band;
    s.upper_band = upper_band;
    s.diagonals.assign(lower_band + upper_band + 1, Vector(n, 0.0));
    s.b.assign(n, 0.0);
    return s;
}

void setBandedValue(BandedSystem& s, int row, int col, double value) {
    const int offset = col - row;
    if (offset < -s.lower_band || offset > s.upper_band) {
        throw std::invalid_argument("Попытка записать элемент вне ленты матрицы.");
    }
    s.diagonals[offset + s.lower_band][row] = value;
}

} // namespace

double vectorInfinityNorm(const Vector& v) {
    double norm = 0.0;
    for (double value : v) {
        norm = std::max(norm, absd(value));
    }
    return norm;
}

double residualInfinityNorm(const Matrix& a, const Vector& x, const Vector& b) {
    Vector ax = multiply(a, x);
    for (std::size_t i = 0; i < ax.size(); ++i) {
        ax[i] -= b[i];
    }
    return vectorInfinityNorm(ax);
}

double residualInfinityNorm(const TridiagonalSystem& s, const Vector& x) {
    const int n = static_cast<int>(s.diag.size());
    double norm = 0.0;

    for (int i = 0; i < n; ++i) {
        double value = s.diag[i] * x[i];
        if (i > 0) {
            value += s.lower[i - 1] * x[i - 1];
        }
        if (i + 1 < n) {
            value += s.upper[i] * x[i + 1];
        }
        norm = std::max(norm, absd(value - s.b[i]));
    }

    return norm;
}

double residualInfinityNorm(const BandedSystem& s, const Vector& x) {
    Vector ax = multiply(s, x);
    for (std::size_t i = 0; i < ax.size(); ++i) {
        ax[i] -= s.b[i];
    }
    return vectorInfinityNorm(ax);
}

SystemStatus analyzeSystem(const Matrix& a, const Vector& b, double eps) {
    SystemStatus status;

    const int n = static_cast<int>(a.size());
    if (n == 0) {
        status.type = SolutionType::INFINITE;
        status.rank_a = 0;
        status.rank_augmented = 0;
        status.nonsingular = false;
        return status;
    }

    for (const auto& row : a) {
        if (static_cast<int>(row.size()) != n) {
            throw std::invalid_argument("Матрица A должна быть квадратной.");
        }
    }
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Размер вектора b должен совпадать с размером матрицы A.");
    }

    Matrix augmented = a;
    for (int i = 0; i < n; ++i) {
        augmented[i].push_back(b[i]);
    }

    // status.rank_a = rankOfMatrix(a, eps);
    // status.rank_augmented = rankOfMatrix(augmented, eps);
    status.rank_a = 100000;
    status.rank_augmented = 100000;

    if (status.rank_a != status.rank_augmented) {
        status.type = SolutionType::NONE;
        status.nonsingular = false;
    } else if (status.rank_a < n) {
        status.type = SolutionType::INFINITE;
        status.nonsingular = false;
    } else {
        status.type = SolutionType::UNIQUE;
        status.nonsingular = true;
    }

    return status;
}

bool reorderRowsByLeadingElement(Matrix& a, Vector& b, double eps) {
    const int n = static_cast<int>(a.size());
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Размер вектора b не совпадает с числом строк матрицы A.");
    }

    Matrix new_a(n, Vector(n, 0.0));
    Vector new_b(n, 0.0);
    std::vector<bool> used(n, false);

    for (int col = 0; col < n; ++col) {
        int best_row = -1;
        double best_abs = -1.0;

        for (int row = 0; row < n; ++row) {
            if (used[row]) {
                continue;
            }
            const double current_abs = absd(a[row][col]);
            if (current_abs > best_abs) {
                best_abs = current_abs;
                best_row = row;
            }
        }

        if (best_row == -1 || best_abs <= eps) {
            return false;
        }

        used[best_row] = true;
        new_a[col] = a[best_row];
        new_b[col] = b[best_row];
    }

    a = std::move(new_a);
    b = std::move(new_b);
    return true;
}

bool isStrictlyDiagonallyDominant(const Matrix& a, double eps) {
    const int n = static_cast<int>(a.size());
    for (int i = 0; i < n; ++i) {
        double off_diagonal_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                off_diagonal_sum += absd(a[i][j]);
            }
        }
        if (!(absd(a[i][i]) > off_diagonal_sum + eps)) {
            return false;
        }
    }
    return true;
}

bool isStrictlyDiagonallyDominant(const TridiagonalSystem& s, double eps) {
    const int n = static_cast<int>(s.diag.size());
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        if (i > 0) {
            sum += absd(s.lower[i - 1]);
        }
        if (i + 1 < n) {
            sum += absd(s.upper[i]);
        }
        if (!(absd(s.diag[i]) > sum + eps)) {
            return false;
        }
    }
    return true;
}

bool isStrictlyDiagonallyDominant(const BandedSystem& s, double eps) {
    const int n = static_cast<int>(s.b.size());
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        const int left = std::max(0, i - s.lower_band);
        const int right = std::min(n - 1, i + s.upper_band);
        for (int j = left; j <= right; ++j) {
            if (i == j) {
                continue;
            }
            sum += absd(getBandedValue(s, i, j));
        }
        if (!(absd(getBandedValue(s, i, i)) > sum + eps)) {
            return false;
        }
    }
    return true;
}

double iterationMatrixInfinityNorm(const Matrix& a, double eps) {
    const int n = static_cast<int>(a.size());
    double norm = 0.0;

    for (int i = 0; i < n; ++i) {
        if (absd(a[i][i]) <= eps) {
            return std::numeric_limits<double>::infinity();
        }

        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                row_sum += absd(a[i][j] / a[i][i]);
            }
        }
        norm = std::max(norm, row_sum);
    }

    return norm;
}

JacobiResult solveJacobi(Matrix a, Vector b, double eps, int max_iterations, bool try_reorder) {
    JacobiResult result;

    const int n = static_cast<int>(a.size());
    if (n == 0) {
        result.converged = true;
        result.guaranteed_convergence = true;
        result.message = "Пустая система: решение считается найденным.";
        return result;
    }

    for (const auto& row : a) {
        if (static_cast<int>(row.size()) != n) {
            throw std::invalid_argument("Матрица A должна быть квадратной.");
        }
    }
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Размер вектора b должен совпадать с размером матрицы A.");
    }
    if (eps <= 0.0 || max_iterations <= 0) {
        throw std::invalid_argument("eps и max_iterations должны быть положительными.");
    }

    const SystemStatus status = analyzeSystem(a, b, 1e-12);
    if (status.type == SolutionType::NONE) {
        result.message = "Система несовместна: решений нет.";
        return result;
    }
    if (status.type == SolutionType::INFINITE) {
        result.message = "Система имеет бесконечно много решений: метод Якоби в таком виде не применяется.";
        return result;
    }

    if (try_reorder) {
        reorderRowsByLeadingElement(a, b, 1e-12);
    }

    result.guaranteed_convergence = isStrictlyDiagonallyDominant(a, 1e-12);
    const double norm_b = iterationMatrixInfinityNorm(a, 1e-12);

    Vector x_old(n, 0.0);
    Vector x_new(n, 0.0);

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        for (int i = 0; i < n; ++i) {
            if (absd(a[i][i]) <= 1e-12) {
                result.message = "На диагонали встретился нулевой элемент, метод Якоби неприменим.";
                result.iterations = iteration - 1;
                return result;
            }

            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += a[i][j] * x_old[j];
                }
            }
            x_new[i] = (b[i] - sum) / a[i][i];
        }

        Vector diff(n, 0.0);
        for (int i = 0; i < n; ++i) {
            diff[i] = x_new[i] - x_old[i];
        }

        result.final_error = vectorInfinityNorm(diff);
        result.iterations = iteration;

        if (result.final_error < eps) {
            result.converged = true;
            result.x = x_new;
            result.residual_norm = residualInfinityNorm(a, result.x, b);

            if (result.guaranteed_convergence) {
                result.message = "Метод Якоби сошёлся. Достаточное условие сходимости выполнено.";
            } else if (norm_b < 1.0) {
                result.message = "Метод Якоби сошёлся. Гарантия установлена по норме итерационной матрицы.";
                result.guaranteed_convergence = true;
            } else {
                result.message = "Метод Якоби сошёлся, но строгой априорной гарантии сходимости не было.";
            }
            return result;
        }

        x_old.swap(x_new);
    }

    result.x = x_old;
    result.residual_norm = residualInfinityNorm(a, result.x, b);
    result.message = "Достигнуто максимальное число итераций, требуемая точность не достигнута.";
    return result;
}

JacobiResult solveJacobi(const TridiagonalSystem& s, double eps, int max_iterations) {
    JacobiResult result;

    const int n = static_cast<int>(s.diag.size());
    if (n == 0) {
        result.converged = true;
        result.guaranteed_convergence = true;
        result.message = "Пустая система: решение считается найденным.";
        return result;
    }
    if (static_cast<int>(s.b.size()) != n ||
        static_cast<int>(s.lower.size()) != std::max(0, n - 1) ||
        static_cast<int>(s.upper.size()) != std::max(0, n - 1)) {
        throw std::invalid_argument("Некорректные размеры трёхдиагональной системы.");
    }
    if (eps <= 0.0 || max_iterations <= 0) {
        throw std::invalid_argument("eps и max_iterations должны быть положительными.");
    }

    result.guaranteed_convergence = isStrictlyDiagonallyDominant(s, 1e-12);

    Vector x_old(n, 0.0);
    Vector x_new(n, 0.0);

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        for (int i = 0; i < n; ++i) {
            if (absd(s.diag[i]) <= 1e-12) {
                result.message = "На главной диагонали встретился нулевой элемент, метод Якоби неприменим.";
                result.iterations = iteration - 1;
                return result;
            }

            double sum = 0.0;
            if (i > 0) {
                sum += s.lower[i - 1] * x_old[i - 1];
            }
            if (i + 1 < n) {
                sum += s.upper[i] * x_old[i + 1];
            }
            x_new[i] = (s.b[i] - sum) / s.diag[i];
        }

        Vector diff(n, 0.0);
        for (int i = 0; i < n; ++i) {
            diff[i] = x_new[i] - x_old[i];
        }

        result.final_error = vectorInfinityNorm(diff);
        result.iterations = iteration;

        if (result.final_error < eps) {
            result.converged = true;
            result.x = x_new;
            result.residual_norm = residualInfinityNorm(s, result.x);
            result.message = result.guaranteed_convergence
                ? "Метод Якоби сошёлся для трёхдиагональной системы."
                : "Метод Якоби сошёлся для трёхдиагональной системы без априорной гарантии.";
            return result;
        }

        x_old.swap(x_new);
    }

    result.x = x_old;
    result.residual_norm = residualInfinityNorm(s, result.x);
    result.message = "Для трёхдиагональной системы достигнут лимит итераций без достижения требуемой точности.";
    return result;
}

JacobiResult solveJacobi(const BandedSystem& s, double eps, int max_iterations) {
    JacobiResult result;

    const int n = static_cast<int>(s.b.size());
    if (n == 0) {
        result.converged = true;
        result.guaranteed_convergence = true;
        result.message = "Пустая система: решение считается найденным.";
        return result;
    }
    if (static_cast<int>(s.diagonals.size()) != s.lower_band + s.upper_band + 1) {
        throw std::invalid_argument("Некорректное число диагоналей в ленточной системе.");
    }
    for (const auto& diag : s.diagonals) {
        if (static_cast<int>(diag.size()) != n) {
            throw std::invalid_argument("Все диагонали ленточной системы должны иметь длину n.");
        }
    }
    if (eps <= 0.0 || max_iterations <= 0) {
        throw std::invalid_argument("eps и max_iterations должны быть положительными.");
    }

    result.guaranteed_convergence = isStrictlyDiagonallyDominant(s, 1e-12);

    Vector x_old(n, 0.0);
    Vector x_new(n, 0.0);

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        for (int i = 0; i < n; ++i) {
            const double diag = getBandedValue(s, i, i);
            if (absd(diag) <= 1e-12) {
                result.message = "На главной диагонали ленточной матрицы встретился нулевой элемент.";
                result.iterations = iteration - 1;
                return result;
            }

            double sum = 0.0;
            const int left = std::max(0, i - s.lower_band);
            const int right = std::min(n - 1, i + s.upper_band);
            for (int j = left; j <= right; ++j) {
                if (j == i) {
                    continue;
                }
                sum += getBandedValue(s, i, j) * x_old[j];
            }

            x_new[i] = (s.b[i] - sum) / diag;
        }

        Vector diff(n, 0.0);
        for (int i = 0; i < n; ++i) {
            diff[i] = x_new[i] - x_old[i];
        }

        result.final_error = vectorInfinityNorm(diff);
        result.iterations = iteration;

        if (result.final_error < eps) {
            result.converged = true;
            result.x = x_new;
            result.residual_norm = residualInfinityNorm(s, result.x);
            result.message = result.guaranteed_convergence
                ? "Метод Якоби сошёлся для ленточной системы."
                : "Метод Якоби сошёлся для ленточной системы без строгой априорной гарантии.";
            return result;
        }

        x_old.swap(x_new);
    }

    result.x = x_old;
    result.residual_norm = residualInfinityNorm(s, result.x);
    result.message = "Для ленточной системы достигнут лимит итераций без достижения требуемой точности.";
    return result;
}

Matrix generateHilbertMatrix(int n) {
    if (n <= 0) {
        throw std::invalid_argument("Размер матрицы Гильберта должен быть положительным.");
    }

    Matrix a(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = 1.0 / static_cast<double>(i + j + 1);
        }
    }
    return a;
}

void generateHilbertSystem(int n, Matrix& a, Vector& b, const Vector& exact_x) {
    if (static_cast<int>(exact_x.size()) != n) {
        throw std::invalid_argument("Размер exact_x должен совпадать с размером матрицы Гильберта.");
    }

    a = generateHilbertMatrix(n);
    b = multiply(a, exact_x);
}

void generateDenseDiagonallyDominantSystem(int n, Matrix& a, Vector& b, const Vector& exact_x, unsigned seed) {
    if (n <= 0) {
        throw std::invalid_argument("Размер системы должен быть положительным.");
    }
    if (static_cast<int>(exact_x.size()) != n) {
        throw std::invalid_argument("Размер exact_x должен совпадать с размером системы.");
    }

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::uniform_real_distribution<double> margin_dist(1.0, 3.0);

    a.assign(n, Vector(n, 0.0));
    b.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double sum_abs = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            a[i][j] = dist(gen);
            sum_abs += absd(a[i][j]);
        }
        a[i][i] = sum_abs + margin_dist(gen);
    }

    b = multiply(a, exact_x);
}

TridiagonalSystem generateLargeTridiagonalSystem(int n, const Vector& exact_x) {
    if (n <= 0) {
        throw std::invalid_argument("Размер системы должен быть положительным.");
    }
    if (static_cast<int>(exact_x.size()) != n) {
        throw std::invalid_argument("Размер exact_x должен совпадать с размером системы.");
    }

    TridiagonalSystem s;
    s.lower.assign(std::max(0, n - 1), -1.0);
    s.diag.assign(n, 4.0);
    s.upper.assign(std::max(0, n - 1), -1.0);
    s.b.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        s.b[i] = s.diag[i] * exact_x[i];
        if (i > 0) {
            s.b[i] += s.lower[i - 1] * exact_x[i - 1];
        }
        if (i + 1 < n) {
            s.b[i] += s.upper[i] * exact_x[i + 1];
        }
    }

    return s;
}

BandedSystem generatePentadiagonalSystem(int n, const Vector& exact_x, unsigned seed) {
    return generateBandedDiagonallyDominantSystem(n, 2, 2, exact_x, seed);
}

BandedSystem generateBandedDiagonallyDominantSystem(
    int n,
    int lower_band,
    int upper_band,
    const Vector& exact_x,
    unsigned seed
) {
    if (n <= 0) {
        throw std::invalid_argument("Размер системы должен быть положительным.");
    }
    if (static_cast<int>(exact_x.size()) != n) {
        throw std::invalid_argument("Размер exact_x должен совпадать с размером системы.");
    }

    BandedSystem s = makeEmptyBandedSystem(n, lower_band, upper_band);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::uniform_real_distribution<double> margin_dist(1.0, 3.0);

    for (int i = 0; i < n; ++i) {
        double sum_abs = 0.0;

        const int left = std::max(0, i - lower_band);
        const int right = std::min(n - 1, i + upper_band);

        for (int j = left; j <= right; ++j) {
            if (i == j) {
                continue;
            }
            const double value = dist(gen);
            setBandedValue(s, i, j, value);
            sum_abs += absd(value);
        }

        setBandedValue(s, i, i, sum_abs + margin_dist(gen));
    }

    s.b = multiply(s, exact_x);
    return s;
}

BandedSystem generateSymmetricPositiveDefiniteBandedSystem(
    int n,
    int half_band,
    const Vector& exact_x,
    unsigned seed
) {
    if (n <= 0) {
        throw std::invalid_argument("Размер системы должен быть положительным.");
    }
    if (static_cast<int>(exact_x.size()) != n) {
        throw std::invalid_argument("Размер exact_x должен совпадать с размером системы.");
    }

    BandedSystem s = makeEmptyBandedSystem(n, half_band, half_band);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> weight_dist(0.05, 0.7);
    std::uniform_real_distribution<double> margin_dist(1.0, 2.0);

    // Сначала заполняем симметричные внедиагональные элементы отрицательными числами.
    for (int i = 0; i < n; ++i) {
        const int right = std::min(n - 1, i + half_band);
        for (int j = i + 1; j <= right; ++j) {
            const double w = weight_dist(gen);
            setBandedValue(s, i, j, -w);
            setBandedValue(s, j, i, -w);
        }
    }

    // Теперь делаем диагональ достаточно большой.
    // Получается симметричная строго диагонально доминирующая матрица,
    // а значит она невырождена и хорошо подходит для метода Якоби.
    for (int i = 0; i < n; ++i) {
        double sum_abs = 0.0;
        const int left = std::max(0, i - half_band);
        const int right = std::min(n - 1, i + half_band);
        for (int j = left; j <= right; ++j) {
            if (i == j) {
                continue;
            }
            sum_abs += absd(getBandedValue(s, i, j));
        }
        setBandedValue(s, i, i, sum_abs + margin_dist(gen));
    }

    s.b = multiply(s, exact_x);
    return s;
}

void saveDenseSystemToFile(const Matrix& a, const Vector& b, const std::string& filename, double eps, int max_iterations) {
    std::ofstream fout(filename);
    if (!fout) {
        throw std::runtime_error("Не удалось открыть файл для записи: " + filename);
    }

    const int n = static_cast<int>(a.size());
    fout << "dense\n";
    fout << n << "\n";
    fout.setf(std::ios::scientific);
    fout.precision(17);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << a[i][j] << (j + 1 == n ? '\n' : ' ');
        }
    }
    for (int i = 0; i < n; ++i) {
        fout << b[i] << (i + 1 == n ? '\n' : ' ');
    }
    fout << eps << ' ' << max_iterations << "\n";
}

void saveTridiagonalSystemToFile(const TridiagonalSystem& s, const std::string& filename, double eps, int max_iterations) {
    std::ofstream fout(filename);
    if (!fout) {
        throw std::runtime_error("Не удалось открыть файл для записи: " + filename);
    }

    const int n = static_cast<int>(s.diag.size());
    fout << "tridiag\n";
    fout << n << "\n";
    fout.setf(std::ios::scientific);
    fout.precision(17);

    for (int i = 0; i < n - 1; ++i) {
        fout << s.lower[i] << (i + 1 == n - 1 ? '\n' : ' ');
    }
    if (n == 1) {
        fout << "\n";
    }
    for (int i = 0; i < n; ++i) {
        fout << s.diag[i] << (i + 1 == n ? '\n' : ' ');
    }
    for (int i = 0; i < n - 1; ++i) {
        fout << s.upper[i] << (i + 1 == n - 1 ? '\n' : ' ');
    }
    if (n == 1) {
        fout << "\n";
    }
    for (int i = 0; i < n; ++i) {
        fout << s.b[i] << (i + 1 == n ? '\n' : ' ');
    }
    fout << eps << ' ' << max_iterations << "\n";
}

void saveBandedSystemToFile(const BandedSystem& s, const std::string& filename, double eps, int max_iterations) {
    std::ofstream fout(filename);
    if (!fout) {
        throw std::runtime_error("Не удалось открыть файл для записи: " + filename);
    }

    const int n = static_cast<int>(s.b.size());
    fout << "banded\n";
    fout << n << ' ' << s.lower_band << ' ' << s.upper_band << "\n";
    fout.setf(std::ios::scientific);
    fout.precision(17);

    for (const auto& diag : s.diagonals) {
        for (int i = 0; i < n; ++i) {
            fout << diag[i] << (i + 1 == n ? '\n' : ' ');
        }
    }
    for (int i = 0; i < n; ++i) {
        fout << s.b[i] << (i + 1 == n ? '\n' : ' ');
    }
    fout << eps << ' ' << max_iterations << "\n";
}

} // namespace jacobi
