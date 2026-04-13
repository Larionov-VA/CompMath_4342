#include "methods.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace {

bool IsZeroRow(const DenseMatrix& matrix,
               std::size_t row,
               double zeroTolerance) {
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        if (std::abs(matrix(row, col)) > zeroTolerance) {
            return false;
        }
    }
    return true;
}

bool HasInconsistentRow(const DenseMatrix& matrix,
                        const std::vector<double>& b,
                        std::size_t startRow,
                        double zeroTolerance) {
    for (std::size_t row = startRow; row < matrix.rows(); ++row) {
        if (IsZeroRow(matrix, row, zeroTolerance) &&
            std::abs(b[row]) > zeroTolerance) {
            return true;
        }
    }
    return false;
}

}  // namespace

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

void DenseMatrix::swapRows(std::size_t first, std::size_t second) {
    if (first == second) {
        return;
    }
    for (std::size_t col = 0; col < cols_; ++col) {
        std::swap((*this)(first, col), (*this)(second, col));
    }
}

GaussianResult SolveGaussian(DenseMatrix A,
                             std::vector<double> b,
                             double zeroTolerance) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }
    if (b.size() != A.rows()) {
        throw std::invalid_argument("Vector size must match matrix size.");
    }

    const std::size_t n = A.rows();
    std::vector<std::size_t> pivotColumns;
    std::size_t pivotRow = 0;
    std::size_t rowSwaps = 0;

    for (std::size_t col = 0; col < n && pivotRow < n; ++col) {
        std::size_t bestRow = pivotRow;
        double bestValue = 0.0;

        for (std::size_t row = pivotRow; row < n; ++row) {
            const double candidate = std::abs(A(row, col));
            if (candidate > bestValue) {
                bestValue = candidate;
                bestRow = row;
            }
        }

        if (bestValue <= zeroTolerance) {
            if (HasInconsistentRow(A, b, pivotRow, zeroTolerance)) {
                GaussianResult result;
                result.status = SolutionStatus::NoSolution;
                result.rank = pivotColumns.size();
                result.rowSwaps = rowSwaps;
                return result;
            }
            continue;
        }

        if (bestRow != pivotRow) {
            A.swapRows(bestRow, pivotRow);
            std::swap(b[bestRow], b[pivotRow]);
            ++rowSwaps;
        }

        const double pivot = A(pivotRow, col);
        for (std::size_t row = pivotRow + 1; row < n; ++row) {
            const double factor = A(row, col) / pivot;
            if (std::abs(factor) <= zeroTolerance) {
                A(row, col) = 0.0;
                continue;
            }

            A(row, col) = 0.0;
            for (std::size_t inner = col + 1; inner < n; ++inner) {
                A(row, inner) -= factor * A(pivotRow, inner);
                if (std::abs(A(row, inner)) <= zeroTolerance) {
                    A(row, inner) = 0.0;
                }
            }
            b[row] -= factor * b[pivotRow];
            if (std::abs(b[row]) <= zeroTolerance) {
                b[row] = 0.0;
            }

            if (IsZeroRow(A, row, zeroTolerance) &&
                std::abs(b[row]) > zeroTolerance) {
                GaussianResult result;
                result.status = SolutionStatus::NoSolution;
                result.rank = pivotColumns.size() + 1;
                result.rowSwaps = rowSwaps;
                return result;
            }
        }

        pivotColumns.push_back(col);
        ++pivotRow;
    }

    if (HasInconsistentRow(A, b, pivotRow, zeroTolerance)) {
        GaussianResult result;
        result.status = SolutionStatus::NoSolution;
        result.rank = pivotColumns.size();
        result.rowSwaps = rowSwaps;
        return result;
    }

    std::vector<double> x(n, 0.0);
    for (std::size_t index = pivotColumns.size(); index > 0; --index) {
        const std::size_t row = index - 1;
        const std::size_t col = pivotColumns[row];

        long double rhs = b[row];
        for (std::size_t inner = col + 1; inner < n; ++inner) {
            rhs -= static_cast<long double>(A(row, inner)) * x[inner];
        }

        x[col] = static_cast<double>(rhs / A(row, col));
    }

    GaussianResult result;
    result.status = (pivotColumns.size() == n)
        ? SolutionStatus::Unique
        : SolutionStatus::InfiniteSolutions;
    result.pivotStrategy = PivotStrategy::PartialByColumn;
    result.x = std::move(x);
    result.rank = pivotColumns.size();
    result.rowSwaps = rowSwaps;
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

namespace {

double ComputeMatrixInfinityNorm(const DenseMatrix& matrix) {
    double maxRowSum = 0.0;

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        long double rowSum = 0.0;
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            rowSum += std::abs(matrix(row, col));
        }
        maxRowSum = std::max(maxRowSum, static_cast<double>(rowSum));
    }

    return maxRowSum;
}

}  // namespace

bool ComputeConditionNumberInf(const DenseMatrix& matrix,
                               double zeroTolerance,
                               double& conditionNumber) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }

    const std::size_t n = matrix.rows();
    const double normA = ComputeMatrixInfinityNorm(matrix);

    std::vector<long double> inverseRowSums(n, 0.0L);

    for (std::size_t col = 0; col < n; ++col) {
        std::vector<double> unitVector(n, 0.0);
        unitVector[col] = 1.0;

        const GaussianResult result =
            SolveGaussian(matrix, unitVector, zeroTolerance);

        if (result.status != SolutionStatus::Unique) {
            return false;
        }

        for (std::size_t row = 0; row < n; ++row) {
            inverseRowSums[row] += std::abs(result.x[row]);
        }
    }

    long double normInvA = 0.0L;
    for (long double rowSum : inverseRowSums) {
        normInvA = std::max(normInvA, rowSum);
    }

    conditionNumber =
        static_cast<double>(static_cast<long double>(normA) * normInvA);

    return std::isfinite(conditionNumber);
}

std::string ToString(PivotStrategy strategy) {
    switch (strategy) {
        case PivotStrategy::PartialByColumn:
            return "частичный выбор по столбцу";
    }
    return "неизвестная стратегия";
}

std::string ToString(SolutionStatus status) {
    switch (status) {
        case SolutionStatus::Unique:
            return "единственное решение";
        case SolutionStatus::InfiniteSolutions:
            return "бесконечно много решений";
        case SolutionStatus::NoSolution:
            return "решений нет";
    }
    return "неизвестный статус";
}

double EstimateDenseMatrixMiB(std::size_t n) {
    const long double bytes =
        static_cast<long double>(n) * n * sizeof(double);
    return static_cast<double>(bytes / (1024.0L * 1024.0L));
}

double EstimateAugmentedMatrixMiB(std::size_t n) {
    const long double bytes =
        static_cast<long double>(n) * (n + 1) * sizeof(double);
    return static_cast<double>(bytes / (1024.0L * 1024.0L));
}

long double EstimateGaussianFlops(std::size_t n) {
    const long double size = static_cast<long double>(n);
    return (2.0L / 3.0L) * size * size * size;
}
