#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include <string>
#include <vector>

namespace jacobi {

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Тип решения системы.
enum class SolutionType {
    UNIQUE,     // Единственное решение.
    NONE,       // Решений нет.
    INFINITE    // Решений бесконечно много.
};

// Информация о системе: ранги, тип решения и признак невырожденности.
struct SystemStatus {
    SolutionType type = SolutionType::NONE;
    int rank_a = 0;
    int rank_augmented = 0;
    bool nonsingular = false;
};

// Результат работы метода Якоби.
struct JacobiResult {
    bool converged = false;                // Сошёлся ли метод.
    bool guaranteed_convergence = false;   // Есть ли достаточное условие сходимости.
    int iterations = 0;                    // Число итераций.
    double final_error = 0.0;              // ||x(k+1) - x(k)||_inf.
    double residual_norm = 0.0;            // ||Ax - b||_inf.
    Vector x;                              // Найденное приближение.
    std::string message;                   // Пояснение для пользователя.
};

// Специальное хранение большой трёхдиагональной системы.
// Такой формат очень экономичен по памяти.
struct TridiagonalSystem {
    Vector lower; // Поддиагональ, размер n-1.
    Vector diag;  // Главная диагональ, размер n.
    Vector upper; // Наддиагональ, размер n-1.
    Vector b;     // Правая часть.
};

// Универсальное хранение ленточной матрицы.
// Матрица хранится не полностью, а только по диагоналям.
// diagonals[index][i] = a[i][i + offset], где
// offset = index - lower_band.
struct BandedSystem {
    int lower_band = 0;              // Сколько диагоналей ниже главной хранится.
    int upper_band = 0;              // Сколько диагоналей выше главной хранится.
    std::vector<Vector> diagonals;   // Всего lower_band + upper_band + 1 диагоналей.
    Vector b;                        // Правая часть.
};

// ---------------- Базовые сервисные функции ----------------

double vectorInfinityNorm(const Vector& v);
double residualInfinityNorm(const Matrix& a, const Vector& x, const Vector& b);
double residualInfinityNorm(const TridiagonalSystem& s, const Vector& x);
double residualInfinityNorm(const BandedSystem& s, const Vector& x);

// ---------------- Анализ системы ----------------

SystemStatus analyzeSystem(const Matrix& a, const Vector& b, double eps = 1e-12);
bool reorderRowsByLeadingElement(Matrix& a, Vector& b, double eps = 1e-12);
bool isStrictlyDiagonallyDominant(const Matrix& a, double eps = 1e-12);
bool isStrictlyDiagonallyDominant(const TridiagonalSystem& s, double eps = 1e-12);
bool isStrictlyDiagonallyDominant(const BandedSystem& s, double eps = 1e-12);
double iterationMatrixInfinityNorm(const Matrix& a, double eps = 1e-12);

// ---------------- Решение методом Якоби ----------------

JacobiResult solveJacobi(
    Matrix a,
    Vector b,
    double eps = 1e-10,
    int max_iterations = 10000,
    bool try_reorder = true
);

JacobiResult solveJacobi(
    const TridiagonalSystem& s,
    double eps = 1e-10,
    int max_iterations = 100000
);

JacobiResult solveJacobi(
    const BandedSystem& s,
    double eps = 1e-10,
    int max_iterations = 100000
);

// ---------------- Генераторы тестовых систем ----------------

Matrix generateHilbertMatrix(int n);
void generateHilbertSystem(int n, Matrix& a, Vector& b, const Vector& exact_x);

// Плотная строго диагонально доминирующая матрица.
void generateDenseDiagonallyDominantSystem(
    int n,
    Matrix& a,
    Vector& b,
    const Vector& exact_x,
    unsigned seed = 123456789u
);

// Большая трёхдиагональная система.
TridiagonalSystem generateLargeTridiagonalSystem(int n, const Vector& exact_x);

// Большая пятидиагональная система.
BandedSystem generatePentadiagonalSystem(
    int n,
    const Vector& exact_x,
    unsigned seed = 123456789u
);

// Общая большая ленточная система.
BandedSystem generateBandedDiagonallyDominantSystem(
    int n,
    int lower_band,
    int upper_band,
    const Vector& exact_x,
    unsigned seed = 123456789u
);

// Симметричная положительно определённая ленточная система.
BandedSystem generateSymmetricPositiveDefiniteBandedSystem(
    int n,
    int half_band,
    const Vector& exact_x,
    unsigned seed = 123456789u
);

// ---------------- Работа с файлами ----------------

void saveDenseSystemToFile(
    const Matrix& a,
    const Vector& b,
    const std::string& filename,
    double eps,
    int max_iterations
);

void saveTridiagonalSystemToFile(
    const TridiagonalSystem& s,
    const std::string& filename,
    double eps,
    int max_iterations
);

void saveBandedSystemToFile(
    const BandedSystem& s,
    const std::string& filename,
    double eps,
    int max_iterations
);

} // namespace jacobi

#endif
