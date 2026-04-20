#include "methods.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

DenseMatrix::DenseMatrix() : rows_(0), cols_(0) {
}

DenseMatrix::DenseMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows),
      cols_(cols),
      data_(rows * cols, 0.0) {
}

std::size_t DenseMatrix::rows() const {
    return rows_;
}

std::size_t DenseMatrix::cols() const {
    return cols_;
}

double& DenseMatrix::operator()(std::size_t row, std::size_t col) {
    return data_[row * cols_ + col];
}

double DenseMatrix::operator()(std::size_t row, std::size_t col) const {
    return data_[row * cols_ + col];
}

JacobiResult SolveJacobi(const DenseMatrix& matrix,
                         const std::vector<double>& b,
                         double eps,
                         std::size_t maxIterations,
                         double zeroTolerance,
                         double divergenceThreshold) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }
    if (b.size() != matrix.rows()) {
        throw std::invalid_argument("Vector size must match matrix size.");
    }

    const std::size_t n = matrix.rows();
    for (std::size_t i = 0; i < n; ++i) {
        if (std::abs(matrix(i, i)) <= zeroTolerance) {
            JacobiResult result;
            result.status = JacobiStatus::ZeroOnDiagonal;
            return result;
        }
    }

    std::vector<double> current(n, 0.0);
    std::vector<double> next(n, 0.0);

    for (std::size_t iteration = 1; iteration <= maxIterations; ++iteration) {
        double stepLinf = 0.0;

        for (std::size_t row = 0; row < n; ++row) {
            long double sum = 0.0;
            for (std::size_t col = 0; col < n; ++col) {
                if (col == row) {
                    continue;
                }
                sum += static_cast<long double>(matrix(row, col)) * current[col];
            }

            next[row] = static_cast<double>(
                (static_cast<long double>(b[row]) - sum) / matrix(row, row)
            );
            stepLinf = std::max(stepLinf, std::abs(next[row] - current[row]));
        }

        if (!std::isfinite(stepLinf) || stepLinf > divergenceThreshold) {
            JacobiResult result;
            result.status = JacobiStatus::Diverged;
            result.x = std::move(next);
            result.iterations = iteration;
            result.lastStepLinf = stepLinf;
            return result;
        }

        current.swap(next);

        if (stepLinf <= eps) {
            JacobiResult result;
            result.status = JacobiStatus::Converged;
            result.x = std::move(current);
            result.iterations = iteration;
            result.lastStepLinf = stepLinf;
            return result;
        }
    }

    JacobiResult result;
    result.status = JacobiStatus::MaxIterationsReached;
    result.x = std::move(current);
    result.iterations = maxIterations;
    result.lastStepLinf = 0.0;
    return result;
}

DenseMatrix MakeMatrix(std::size_t rows,
                       std::size_t cols,
                       const std::vector<double>& values) {
    if (rows * cols != values.size()) {
        throw std::invalid_argument("Row-major value count is incorrect.");
    }

    DenseMatrix matrix(rows, cols);
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            matrix(row, col) = values[row * cols + col];
        }
    }
    return matrix;
}

DenseMatrix MakeHilbertMatrix(std::size_t n) {
    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            matrix(row, col) = 1.0 / static_cast<double>(row + col + 1);
        }
    }
    return matrix;
}

DenseMatrix MakeRandomDenseMatrix(std::size_t n,
                                  std::uint64_t seed,
                                  double minValue,
                                  double maxValue) {
    if (!(minValue < maxValue)) {
        throw std::invalid_argument("Random interval must be non-empty.");
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> distribution(minValue, maxValue);

    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            double value = distribution(rng);
            if (std::abs(value) < 0.1) {
                value = (value < 0.0) ? -0.1 : 0.1;
            }
            matrix(row, col) = value;
        }
    }
    return matrix;
}

DenseMatrix MakeRandomNearlySingularMatrix(std::size_t n,
                                           std::uint64_t seed,
                                           double perturbation) {
    if (perturbation <= 0.0) {
        throw std::invalid_argument("Perturbation must be positive.");
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> valueDistribution(-2.0, 2.0);
    std::uniform_real_distribution<double> scaleDistribution(0.5, 2.0);
    std::uniform_real_distribution<double> noiseDistribution(-1.0, 1.0);

    std::vector<double> baseRow(n, 0.0);
    for (double& value : baseRow) {
        value = valueDistribution(rng);
        if (std::abs(value) < 0.2) {
            value = (value < 0.0) ? -0.2 : 0.2;
        }
    }

    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        const double scale = scaleDistribution(rng);
        for (std::size_t col = 0; col < n; ++col) {
            matrix(row, col) =
                scale * baseRow[col] + perturbation * noiseDistribution(rng);
        }
        matrix(row, row) += perturbation * static_cast<double>(row + 1);
    }
    return matrix;
}

DenseMatrix MakeRandomDiagonalDominantMatrix(std::size_t n,
                                             std::uint64_t seed,
                                             double minValue,
                                             double maxValue) {
    if (!(minValue < maxValue)) {
        throw std::invalid_argument("Random interval must be non-empty.");
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> distribution(minValue, maxValue);
    std::uniform_real_distribution<double> diagonalGap(1.0, 3.0);

    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        double rowSum = 0.0;
        for (std::size_t col = 0; col < n; ++col) {
            if (row == col) {
                continue;
            }

            double value = distribution(rng);
            if (std::abs(value) < 0.05) {
                value = (value < 0.0) ? -0.05 : 0.05;
            }

            matrix(row, col) = value;
            rowSum += std::abs(value);
        }

        const double sign = (distribution(rng) >= 0.0) ? 1.0 : -1.0;
        matrix(row, row) = sign * (rowSum + diagonalGap(rng));
    }
    return matrix;
}

DenseMatrix MakeRandomRankDeficientMatrix(std::size_t n,
                                          std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> valueDistribution(-3.0, 3.0);
    std::uniform_real_distribution<double> scaleDistribution(-2.0, 2.0);

    std::vector<double> baseRow(n, 0.0);
    for (double& value : baseRow) {
        value = valueDistribution(rng);
        if (std::abs(value) < 0.2) {
            value = (value < 0.0) ? -0.2 : 0.2;
        }
    }

    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        double scale = scaleDistribution(rng);
        if (std::abs(scale) < 0.3) {
            scale = (scale < 0.0) ? -1.0 : 1.0;
        }
        for (std::size_t col = 0; col < n; ++col) {
            matrix(row, col) = scale * baseRow[col];
        }
    }
    return matrix;
}

DenseMatrix MakeRandomZeroDiagonalMatrix(std::size_t n,
                                         std::uint64_t seed,
                                         double minValue,
                                         double maxValue) {
    DenseMatrix matrix = MakeRandomDenseMatrix(n, seed, minValue, maxValue);
    for (std::size_t i = 0; i < n; ++i) {
        matrix(i, i) = 0.0;
    }
    return matrix;
}

std::vector<double> MakeRandomVector(std::size_t n,
                                     std::uint64_t seed,
                                     double minValue,
                                     double maxValue) {
    if (!(minValue < maxValue)) {
        throw std::invalid_argument("Random interval must be non-empty.");
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> distribution(minValue, maxValue);

    std::vector<double> values(n, 0.0);
    for (double& value : values) {
        value = distribution(rng);
    }
    return values;
}

std::vector<double> Multiply(const DenseMatrix& matrix,
                             const std::vector<double>& x) {
    if (matrix.cols() != x.size()) {
        throw std::invalid_argument("Matrix/vector dimensions do not match.");
    }

    std::vector<double> result(matrix.rows(), 0.0);
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        long double sum = 0.0;
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            sum += static_cast<long double>(matrix(row, col)) * x[col];
        }
        result[row] = static_cast<double>(sum);
    }
    return result;
}

std::vector<double> Subtract(const std::vector<double>& left,
                             const std::vector<double>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("Vector dimensions do not match.");
    }

    std::vector<double> result(left.size(), 0.0);
    for (std::size_t index = 0; index < left.size(); ++index) {
        result[index] = left[index] - right[index];
    }
    return result;
}

VectorNorms ComputeNorms(const std::vector<double>& values) {
    long double l1 = 0.0;
    long double l2Squared = 0.0;
    double linf = 0.0;

    for (double value : values) {
        const double absolute = std::abs(value);
        l1 += absolute;
        l2Squared += static_cast<long double>(absolute) * absolute;
        linf = std::max(linf, absolute);
    }

    VectorNorms norms{};
    norms.l1 = static_cast<double>(l1);
    norms.l2 = std::sqrt(static_cast<double>(l2Squared));
    norms.linf = linf;
    return norms;
}

DominanceInfo AnalyzeDiagonalDominance(const DenseMatrix& matrix,
                                       double zeroTolerance) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }

    DominanceInfo info;
    info.weak = true;
    info.strict = true;
    info.maxRatio = 0.0;

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        const double diagAbs = std::abs(matrix(row, row));
        if (diagAbs <= zeroTolerance) {
            info.zeroDiagonal = true;
            info.weak = false;
            info.strict = false;
            info.maxRatio = std::numeric_limits<double>::infinity();
            return info;
        }

        long double offDiagonalSum = 0.0;
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            if (col == row) {
                continue;
            }
            offDiagonalSum += std::abs(matrix(row, col));
        }

        const double ratio = static_cast<double>(offDiagonalSum / diagAbs);
        info.maxRatio = std::max(info.maxRatio, ratio);

        if (diagAbs + zeroTolerance < static_cast<double>(offDiagonalSum)) {
            info.weak = false;
            info.strict = false;
        } else if (diagAbs <= static_cast<double>(offDiagonalSum) + zeroTolerance) {
            info.strict = false;
        }
    }

    return info;
}

double EstimateDenseMatrixMiB(std::size_t n) {
    const long double bytes =
        static_cast<long double>(n) * n * sizeof(double);
    return static_cast<double>(bytes / (1024.0L * 1024.0L));
}

long double EstimateJacobiIterationFlops(std::size_t n) {
    const long double size = static_cast<long double>(n);
    return 2.0L * size * size;
}

std::string ToString(JacobiStatus status) {
    switch (status) {
        case JacobiStatus::Converged:
            return "сошелся";
        case JacobiStatus::MaxIterationsReached:
            return "достигнут лимит итераций";
        case JacobiStatus::ZeroOnDiagonal:
            return "нулевой элемент на диагонали";
        case JacobiStatus::Diverged:
            return "расходится";
    }
    return "неизвестный статус";
}
