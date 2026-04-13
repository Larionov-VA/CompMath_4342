#include <clocale>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "methods.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {

constexpr double kZeroTolerance = 1e-12;
constexpr std::size_t kVeryLargeN = 10000;
constexpr std::size_t kFullMatrixLimit = 10;
constexpr std::size_t kPreviewRows = 6;
constexpr std::size_t kPreviewCols = 6;
constexpr std::size_t kPreviewVectorItems = 8;
constexpr std::size_t kConditionNumberLimit = 30;

struct ExperimentReport {
    std::string name;
    std::string matrixType;
    std::string note;
    std::string systemPreview;
    std::uint64_t seed = 0;
    std::size_t n = 0;
    bool executed = false;
    SolutionStatus status = SolutionStatus::NoSolution;
    PivotStrategy pivotStrategy = PivotStrategy::PartialByColumn;
    std::size_t rank = 0;
    std::size_t rowSwaps = 0;
    double elapsedMs = 0.0;
    double conditionNumberInf = std::numeric_limits<double>::quiet_NaN();
    std::string conditionNote;
    VectorNorms solutionError{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    VectorNorms residual{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
};

void SetupConsole() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    std::setlocale(LC_ALL, ".UTF-8");
}

bool HasFlag(int argc, char* argv[], const std::string& flag) {
    for (int index = 1; index < argc; ++index) {
        if (flag == argv[index]) {
            return true;
        }
    }
    return false;
}

std::uint64_t GetRunSeed(int argc, char* argv[]) {
    const std::string prefix = "--seed=";
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind(prefix, 0) == 0) {
            return std::stoull(argument.substr(prefix.size()));
        }
        if (argument == "--seed" && index + 1 < argc) {
            return std::stoull(argv[index + 1]);
        }
    }

    std::random_device device;
    const std::uint64_t high =
        static_cast<std::uint64_t>(device()) << 32;
    const std::uint64_t low = static_cast<std::uint64_t>(device());
    return high ^ low;
}

std::uint64_t DeriveSeed(std::uint64_t seed, std::uint64_t salt) {
    return seed + 0x9E3779B97F4A7C15ULL * salt;
}

DenseMatrix MakeRandomPivotRequiredMatrix(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> distribution(-5.0, 5.0);

    for (int attempt = 0; attempt < 2000; ++attempt) {
        DenseMatrix matrix = MakeRandomDenseMatrix(3, rng());
        matrix(0, 0) = 0.0;

        if (std::abs(matrix(1, 0)) < 0.5) {
            matrix(1, 0) = (distribution(rng) >= 0.0) ? 1.0 : -1.0;
        }
        if (std::abs(matrix(2, 0)) < 0.5) {
            matrix(2, 0) = (distribution(rng) >= 0.0) ? 2.0 : -2.0;
        }

        const std::vector<double> xExact = MakeRandomVector(
            3,
            DeriveSeed(seed, static_cast<std::uint64_t>(attempt) + 1)
        );
        const GaussianResult probe = SolveGaussian(
            matrix,
            Multiply(matrix, xExact),
            kZeroTolerance
        );
        if (probe.status == SolutionStatus::Unique) {
            return matrix;
        }
    }

    throw std::runtime_error(
        "Не удалось сгенерировать случайную матрицу с обязательной перестановкой строк."
    );
}

DenseMatrix MakeRandomRankDeficientMatrix(std::size_t n, std::uint64_t seed) {
    const std::vector<double> baseRow = MakeRandomVector(n, DeriveSeed(seed, 1));
    std::mt19937_64 rng(DeriveSeed(seed, 2));
    std::uniform_real_distribution<double> scaleDistribution(-3.0, 3.0);

    DenseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        double scale = scaleDistribution(rng);
        if (std::abs(scale) < 0.5) {
            scale = (scale < 0.0) ? -1.0 : 1.0;
        }
        for (std::size_t col = 0; col < n; ++col) {
            matrix(row, col) = scale * baseRow[col];
        }
    }
    return matrix;
}

std::string FormatMatrix(const DenseMatrix& matrix) {
    const bool fullMatrix = matrix.rows() <= kFullMatrixLimit &&
                            matrix.cols() <= kFullMatrixLimit;
    const std::size_t shownRows = fullMatrix
        ? matrix.rows()
        : std::min(matrix.rows(), kPreviewRows);
    const std::size_t shownCols = fullMatrix
        ? matrix.cols()
        : std::min(matrix.cols(), kPreviewCols);

    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    if (fullMatrix) {
        out << "Матрица A:\n";
    } else {
        out << "Матрица A, показан левый верхний фрагмент "
            << shownRows << "x" << shownCols
            << " из " << matrix.rows() << "x" << matrix.cols() << ":\n";
    }

    for (std::size_t row = 0; row < shownRows; ++row) {
        out << "  [";
        for (std::size_t col = 0; col < shownCols; ++col) {
            out << std::setw(12) << matrix(row, col);
            if (col + 1 < shownCols) {
                out << ' ';
            }
        }
        if (!fullMatrix && shownCols < matrix.cols()) {
            out << " ...";
        }
        out << " ]\n";
    }

    if (!fullMatrix && shownRows < matrix.rows()) {
        out << "  ...\n";
    }

    return out.str();
}

std::string FormatVector(const std::vector<double>& values,
                         const std::string& title) {
    const bool fullVector = values.size() <= kFullMatrixLimit;
    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    if (fullVector) {
        out << title << " = [";
        for (std::size_t index = 0; index < values.size(); ++index) {
            out << values[index];
            if (index + 1 < values.size()) {
                out << ", ";
            }
        }
        out << "]\n";
        return out.str();
    }

    const std::size_t shown = std::min(values.size(), kPreviewVectorItems);
    out << title << ", первые " << shown
        << " элементов из " << values.size() << ": [";
    for (std::size_t index = 0; index < shown; ++index) {
        out << values[index];
        if (index + 1 < shown) {
            out << ", ";
        }
    }
    out << "]\n";
    return out.str();
}

std::string BuildSystemPreview(const DenseMatrix& matrix,
                               const std::vector<double>& b,
                               const std::vector<double>* exactX) {
    std::ostringstream out;
    out << FormatMatrix(matrix);
    out << FormatVector(b, "Вектор b");
    if (exactX != nullptr) {
        out << FormatVector(*exactX, "Точное решение x*");
    }
    return out.str();
}

ExperimentReport RunExperiment(const std::string& name,
                               const std::string& matrixType,
                               DenseMatrix matrix,
                               std::vector<double> b,
                               const std::vector<double>* exactX,
                               const std::string& note,
                               std::uint64_t seed) {
    const DenseMatrix originalMatrix = matrix;
    const std::vector<double> originalB = b;

    const auto start = std::chrono::steady_clock::now();
    const GaussianResult result =
        SolveGaussian(std::move(matrix), std::move(b), kZeroTolerance);
    const auto finish = std::chrono::steady_clock::now();

    ExperimentReport report;
    report.name = name;
    report.matrixType = matrixType;
    report.note = note;
    report.systemPreview = BuildSystemPreview(originalMatrix, originalB, exactX);
    report.seed = seed;
    report.n = originalMatrix.rows();
    report.executed = true;
    report.status = result.status;
    report.pivotStrategy = result.pivotStrategy;
    report.rank = result.rank;
    report.rowSwaps = result.rowSwaps;
    report.elapsedMs =
        std::chrono::duration<double, std::milli>(finish - start).count();

    if (!result.x.empty()) {
        report.residual =
            ComputeNorms(Subtract(Multiply(originalMatrix, result.x), originalB));
    }

    if (exactX != nullptr && result.status == SolutionStatus::Unique) {
        report.solutionError = ComputeNorms(Subtract(result.x, *exactX));
    }

        if (result.status == SolutionStatus::Unique) {
        if (report.n <= kConditionNumberLimit) {
            double conditionNumber = 0.0;
            if (ComputeConditionNumberInf(
                    originalMatrix,
                    kZeroTolerance,
                    conditionNumber)) {
                report.conditionNumberInf = conditionNumber;
            } else {
                report.conditionNote =
                    "не удалось вычислить cond_inf(A)";
            }
        } else {
            report.conditionNote =
                "не вычислялось: для больших матриц расчет слишком дорогой";
        }
    } else {
        report.conditionNote =
            "не вычисляется для вырожденной матрицы";
    }

    return report;
}

ExperimentReport RunExperimentWithExactSolution(const std::string& name,
                                                const std::string& matrixType,
                                                DenseMatrix matrix,
                                                const std::vector<double>& exactX,
                                                const std::string& note,
                                                std::uint64_t seed) {
    const std::vector<double> b = Multiply(matrix, exactX);
    return RunExperiment(name, matrixType, std::move(matrix), b, &exactX, note, seed);
}

ExperimentReport MakeSkippedLargeReport(std::uint64_t seed) {
    ExperimentReport report;
    report.name = "Случайная диагонально доминирующая матрица 10000x10000";
    report.matrixType = "случайная диагонально доминирующая";
    report.note =
        "По умолчанию тест пропущен. Запустите программу с флагом --run-10000, чтобы принудительно выполнить этот эксперимент.";
    report.seed = seed;
    report.n = kVeryLargeN;
    report.executed = false;
    return report;
}

void PrintNormLine(const std::string& label, const VectorNorms& norms) {
    std::cout << "  " << label
              << ": L1 = " << std::scientific << norms.l1
              << ", L2 = " << norms.l2
              << ", Linf = " << norms.linf << '\n';
}

void PrintReport(const ExperimentReport& report) {
    std::cout << "Случай: " << report.name << '\n';
    std::cout << "  Размер n = " << report.n
              << ", тип матрицы = " << report.matrixType << '\n';
    if (report.seed != 0) {
        std::cout << "  Seed = " << report.seed << '\n';
    }
    if (!report.systemPreview.empty()) {
        std::cout << report.systemPreview;
    }

    if (!report.executed) {
        std::cout << "  Статус = пропущен\n";
        std::cout << "  Оценка памяти только под A, МиБ = "
                  << std::fixed << std::setprecision(3)
                  << EstimateDenseMatrixMiB(report.n) << '\n';
        std::cout << "  Оценка памяти под [A|b], МиБ = "
                  << EstimateAugmentedMatrixMiB(report.n) << '\n';
        std::cout << "  Оценка числа операций Гаусса = "
                  << std::scientific << EstimateGaussianFlops(report.n) << '\n';
        std::cout << "  Примечание = " << report.note << "\n\n";
        return;
    }

    std::cout << "  Статус = " << ToString(report.status)
              << ", ранг = " << report.rank << "/" << report.n
              << ", перестановок строк = " << report.rowSwaps
              << ", время, мс = " << std::fixed << std::setprecision(3)
              << report.elapsedMs << '\n';
    std::cout << "  Стратегия выбора ведущего элемента = "
              << ToString(report.pivotStrategy) << '\n';
    PrintNormLine("невязка", report.residual);
        if (!std::isnan(report.conditionNumberInf)) {
        std::cout << "  cond_inf(A) = "
                  << std::scientific << report.conditionNumberInf << '\n';
    } else if (!report.conditionNote.empty()) {
        std::cout << "  cond_inf(A): " << report.conditionNote << '\n';
    }

    if (!std::isnan(report.solutionError.l2)) {
        PrintNormLine("ошибка решения", report.solutionError);
    }

    if (!report.note.empty()) {
        std::cout << "  Примечание = " << report.note << '\n';
    }
    std::cout << '\n';
}

void WriteCsv(const std::string& filename,
              const std::vector<ExperimentReport>& reports) {
    std::ofstream out(filename);
    out << "случай;seed;тип_матрицы;n;выполнен;статус;ранг;перестановки_строк;время_мс;"
       "cond_inf;комментарий_cond;ошибка_l1;ошибка_l2;ошибка_linf;"
       "невязка_l1;невязка_l2;невязка_linf;примечание\n";
    out << std::setprecision(17);

    for (const auto& report : reports) {
        out << report.name << ';'
            << report.seed << ';'
            << report.matrixType << ';'
            << report.n << ';'
            << (report.executed ? 1 : 0) << ';'
            << (report.executed ? ToString(report.status) : "пропущен") << ';'
            << report.rank << ';'
            << report.rowSwaps << ';'
            << report.elapsedMs << ';'
            << report.conditionNumberInf << ';'
            << report.conditionNote << ';'
            << report.solutionError.l1 << ';'
            << report.solutionError.l2 << ';'
            << report.solutionError.linf << ';'
            << report.residual.l1 << ';'
            << report.residual.l2 << ';'
            << report.residual.linf << ';'
            << report.note << '\n';
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        SetupConsole();
        const bool runVeryLarge = HasFlag(argc, argv, "--run-10000");
        const std::uint64_t runSeed = GetRunSeed(argc, argv);
        std::mt19937_64 seedGenerator(runSeed);
        const auto nextExperimentSeed = [&seedGenerator]() {
            return seedGenerator();
        };

        std::vector<ExperimentReport> reports;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Решение СЛАУ методом Гаусса\n";
        std::cout << "Стратегия выбора ведущего элемента: частичный выбор по столбцу\n";
        std::cout << "Порог нуля: " << std::scientific
                  << kZeroTolerance << "\n";
        std::cout << "Базовый seed запуска: " << runSeed << "\n\n";

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomDiagonalDominantMatrix(
                3,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                3,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайный базовый пример 3x3",
                "случайная диагонально доминирующая",
                matrix,
                xExact,
                "Матрица и точное решение сгенерированы случайно.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomPivotRequiredMatrix(
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                3,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайный пример 3x3 с обязательной перестановкой строк",
                "случайная плотная",
                matrix,
                xExact,
                "Первый элемент в столбце занулен специально, поэтому требуется перестановка строк.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomRankDeficientMatrix(
                3,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                3,
                DeriveSeed(seed, 2)
            );
            const std::vector<double> b = Multiply(matrix, xExact);
            reports.push_back(RunExperiment(
                "Случайная вырожденная совместная система 3x3",
                "случайная вырожденная плотная",
                matrix,
                b,
                nullptr,
                "Строки матрицы линейно зависимы, поэтому система совместна, но решение не единственное.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomRankDeficientMatrix(
                3,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                3,
                DeriveSeed(seed, 2)
            );
            std::vector<double> b = Multiply(matrix, xExact);
            b.back() += 1.0;
            reports.push_back(RunExperiment(
                "Случайная вырожденная несовместная система 3x3",
                "случайная вырожденная плотная",
                matrix,
                b,
                nullptr,
                "Правая часть специально искажена, поэтому возникает противоречивая строка.",
                seed
            ));
        }

        for (const std::size_t n : std::vector<std::size_t>{5, 8, 10}) {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomNearlySingularMatrix(
                n,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                n,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная почти вырожденная матрица n=" + std::to_string(n),
                "случайная плохо обусловленная",
                matrix,
                xExact,
                "Матрица построена как почти линейно зависимая, чтобы показать влияние плохой обусловленности.",
                seed
            ));
        }

        for (const std::size_t n : std::vector<std::size_t>{5, 8, 20}) {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeHilbertMatrix(n);
            const std::vector<double> xExact(n, 1.0);

            reports.push_back(RunExperimentWithExactSolution(
                "Матрица Гильберта n=" + std::to_string(n),
                "матрица Гильберта",
                matrix,
                xExact,
                "Классический плохо обусловленный тест, точное решение выбрано как x* = (1, ..., 1).",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();

            std::size_t customN;
            std::cout << "Введите размер матрицы: ";
            std::cin >> customN;

            const DenseMatrix matrix = MakeRandomDenseMatrix(
                customN,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                customN,
                DeriveSeed(seed, 2)
            );

            reports.push_back(RunExperimentWithExactSolution(
                "Пользовательская матрица " + std::to_string(customN) + "x" + std::to_string(customN),
                "случайная плотная",
                matrix,
                xExact,
                "Размер матрицы введен пользователем.",
                seed
            ));
        }

        const std::uint64_t largeSeed = nextExperimentSeed();
        if (runVeryLarge) {
            const DenseMatrix matrix = MakeRandomDiagonalDominantMatrix(
                kVeryLargeN,
                DeriveSeed(largeSeed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                kVeryLargeN,
                DeriveSeed(largeSeed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная диагонально доминирующая матрица 10000x10000",
                "случайная диагонально доминирующая",
                matrix,
                xExact,
                "Эксперимент был явно включен через флаг --run-10000.",
                largeSeed
            ));
        } else {
            reports.push_back(MakeSkippedLargeReport(largeSeed));
        }

        for (const auto& report : reports) {
            PrintReport(report);
        }

        WriteCsv("gauss_results.csv", reports);
        std::cout << "CSV-отчет сохранен в файл gauss_results.csv\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Ошибка: " << error.what() << '\n';
        return 1;
    }
}
