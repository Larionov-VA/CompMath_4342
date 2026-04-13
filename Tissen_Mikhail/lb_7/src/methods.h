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

    void swapRows(std::size_t first, std::size_t second);

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;
};

enum class PivotStrategy {
    PartialByColumn
};

enum class SolutionStatus {
    Unique,
    InfiniteSolutions,
    NoSolution
};

struct VectorNorms {
    double l1;
    double l2;
    double linf;
};

struct GaussianResult {
    SolutionStatus status = SolutionStatus::NoSolution;
    PivotStrategy pivotStrategy = PivotStrategy::PartialByColumn;
    std::vector<double> x;
    std::size_t rank = 0;
    std::size_t rowSwaps = 0;
};

GaussianResult SolveGaussian(DenseMatrix A,
                             std::vector<double> b,
                             double zeroTolerance = 1e-12);

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

std::vector<double> Multiply(const DenseMatrix& matrix,
                             const std::vector<double>& x);
std::vector<double> Subtract(const std::vector<double>& left,
                             const std::vector<double>& right);

std::vector<double> MakeRandomVector(std::size_t n,
                                     std::uint64_t seed,
                                     double minValue = -5.0,
                                     double maxValue = 5.0);

VectorNorms ComputeNorms(const std::vector<double>& values);

bool ComputeConditionNumberInf(const DenseMatrix& matrix,
                               double zeroTolerance,
                               double& conditionNumber);

std::string ToString(PivotStrategy strategy);
std::string ToString(SolutionStatus status);

double EstimateDenseMatrixMiB(std::size_t n);
double EstimateAugmentedMatrixMiB(std::size_t n);
long double EstimateGaussianFlops(std::size_t n);