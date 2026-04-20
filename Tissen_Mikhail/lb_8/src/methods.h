#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class DenseMatrix {
public:
    DenseMatrix();
    DenseMatrix(std::size_t rows, std::size_t cols);

    std::size_t rows() const;
    std::size_t cols() const;

    double& operator()(std::size_t row, std::size_t col);
    double operator()(std::size_t row, std::size_t col) const;

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;
};

enum class JacobiStatus {
    Converged,
    MaxIterationsReached,
    ZeroOnDiagonal,
    Diverged
};

struct VectorNorms {
    double l1;
    double l2;
    double linf;
};

struct DominanceInfo {
    bool weak = false;
    bool strict = false;
    bool zeroDiagonal = false;
    double maxRatio = 0.0;
};

struct JacobiResult {
    JacobiStatus status = JacobiStatus::MaxIterationsReached;
    std::vector<double> x;
    std::size_t iterations = 0;
    double lastStepLinf = 0.0;
};

JacobiResult SolveJacobi(const DenseMatrix& matrix,
                         const std::vector<double>& b,
                         double eps = 1e-10,
                         std::size_t maxIterations = 10000,
                         double zeroTolerance = 1e-12,
                         double divergenceThreshold = 1e100);

DenseMatrix MakeMatrix(std::size_t rows,
                       std::size_t cols,
                       const std::vector<double>& values);

DenseMatrix MakeHilbertMatrix(std::size_t n);
DenseMatrix MakeRandomDenseMatrix(std::size_t n,
                                  std::uint64_t seed,
                                  double minValue = -10.0,
                                  double maxValue = 10.0);
DenseMatrix MakeRandomNearlySingularMatrix(std::size_t n,
                                           std::uint64_t seed,
                                           double perturbation = 1e-6);
DenseMatrix MakeRandomDiagonalDominantMatrix(std::size_t n,
                                             std::uint64_t seed,
                                             double minValue = -3.0,
                                             double maxValue = 3.0);
DenseMatrix MakeRandomRankDeficientMatrix(std::size_t n,
                                          std::uint64_t seed);
DenseMatrix MakeRandomZeroDiagonalMatrix(std::size_t n,
                                         std::uint64_t seed,
                                         double minValue = -5.0,
                                         double maxValue = 5.0);

std::vector<double> MakeRandomVector(std::size_t n,
                                     std::uint64_t seed,
                                     double minValue = -5.0,
                                     double maxValue = 5.0);
std::vector<double> Multiply(const DenseMatrix& matrix,
                             const std::vector<double>& x);
std::vector<double> Subtract(const std::vector<double>& left,
                             const std::vector<double>& right);
VectorNorms ComputeNorms(const std::vector<double>& values);
DominanceInfo AnalyzeDiagonalDominance(const DenseMatrix& matrix,
                                       double zeroTolerance = 1e-12);

double EstimateDenseMatrixMiB(std::size_t n);
long double EstimateJacobiIterationFlops(std::size_t n);

std::string ToString(JacobiStatus status);
