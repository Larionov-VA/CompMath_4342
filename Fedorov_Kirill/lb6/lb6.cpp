    #include <algorithm>
    #include <chrono>
    #include <clocale>
    #include <cstddef>
    #include <cmath>
    #include <iomanip>
    #include <iostream>
    #include <limits>
    #include <random>
    #include <sstream>
    #include <stdexcept>
    #include <string>
    #include <utility>
    #include <vector>

    class SparseMatrix {
    public:
        explicit SparseMatrix(std::size_t n = 0) : rows_(n) {}

        std::size_t size() const { return rows_.size(); }

        void set(std::size_t i, std::size_t j, double v) {
            if (i >= rows_.size() || j >= rows_.size()) {
                throw std::out_of_range("SparseMatrix::set index out of range");
            }
            auto& r = rows_[i];
            auto it = std::find_if(r.begin(), r.end(), [j](const auto& p) { return p.first == j; });
            if (std::abs(v) <= 1e-15) {
                if (it != r.end()) {
                    r.erase(it);
                }
                return;
            }
            if (it != r.end()) {
                it->second = v;
            } else {
                r.emplace_back(j, v);
            }
        }

        double get(std::size_t i, std::size_t j) const {
            if (i >= rows_.size() || j >= rows_.size()) {
                throw std::out_of_range("SparseMatrix::get index out of range");
            }
            const auto& r = rows_[i];
            auto it = std::find_if(r.begin(), r.end(), [j](const auto& p) { return p.first == j; });
            return it == r.end() ? 0.0 : it->second;
        }

        const std::vector<std::pair<std::size_t, double>>& row(std::size_t i) const {
            return rows_[i];
        }

        void swapRows(std::size_t i, std::size_t j) {
            if (i == j) {
                return;
            }
            std::swap(rows_[i], rows_[j]);
        }

        std::size_t memoryBytes() const {
            std::size_t bytes = sizeof(*this);
            bytes += rows_.capacity() * sizeof(std::vector<std::pair<std::size_t, double>>);
            for (const auto& r : rows_) {
                bytes += r.capacity() * sizeof(std::pair<std::size_t, double>);
            }
            return bytes;
        }

    private:
        std::vector<std::vector<std::pair<std::size_t, double>>> rows_;
    };

    enum class SolutionType { UNIQUE, INFINITE, NONE, UNDETERMINED };

    struct JacobiConfig {
        double tolerance = 1e-10;
        std::size_t maxIterations = 100000;
        double zeroEps = 1e-12;
        bool usePivot = true;
        std::size_t rankCheckLimit = 200;
    };

    struct SolveResult {
        bool success = false;
        bool converged = false;
        SolutionType type = SolutionType::UNDETERMINED;
        int rankA = -1;
        int rankAb = -1;
        std::size_t iterations = 0;
        double achievedDelta = 0.0;
        double residualInf = std::numeric_limits<double>::infinity();
        std::size_t memoryBytes = 0;
        std::vector<double> x;
        std::string message;
    };

    struct RankResult {
        SolutionType type = SolutionType::UNDETERMINED;
        int rankA = -1;
        int rankAb = -1;
    };

    static std::vector<double> multiply(const SparseMatrix& a, const std::vector<double>& x) {
        std::vector<double> out(a.size(), 0.0);
        for (std::size_t i = 0; i < a.size(); ++i) {
            for (const auto& [j, v] : a.row(i)) {
                out[i] += v * x[j];
            }
        }
        return out;
    }

    static double maxAbsError(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            return std::numeric_limits<double>::infinity();
        }
        double m = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            m = std::max(m, std::abs(a[i] - b[i]));
        }
        return m;
    }

    static double l2ErrorNorm(const std::vector<double>& xRef, const std::vector<double>& xApprox) {
        if (xRef.size() != xApprox.size()) {
            throw std::invalid_argument("Размеры векторов для нормы ошибки не совпадают");
        }
        long double sum = 0.0L;
        for (std::size_t i = 0; i < xRef.size(); ++i) {
            const long double d = static_cast<long double>(xRef[i]) - static_cast<long double>(xApprox[i]);
            sum += d * d;
        }
        return std::sqrt(static_cast<double>(sum));
    }

    static std::string formatScientific(double value, int precision = 6) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision) << value;
        return oss.str();
    }

    static double residualInfNorm(const SparseMatrix& a, const std::vector<double>& x, const std::vector<double>& b) {
        const std::vector<double> ax = multiply(a, x);
        double m = 0.0;
        for (std::size_t i = 0; i < b.size(); ++i) {
            m = std::max(m, std::abs(ax[i] - b[i]));
        }
        return m;
    }

    static bool isStrictlyDiagonallyDominant(const SparseMatrix& a, double eps) {
        for (std::size_t i = 0; i < a.size(); ++i) {
            double diag = 0.0;
            double sum = 0.0;
            for (const auto& [j, v] : a.row(i)) {
                if (j == i) {
                    diag = std::abs(v);
                } else {
                    sum += std::abs(v);
                }
            }
            if (!(diag > sum + eps)) {
                return false;
            }
        }
        return true;
    }

    static bool hasInconsistentZeroRow(const SparseMatrix& a, const std::vector<double>& b, double eps) {
        for (std::size_t i = 0; i < a.size(); ++i) {
            bool zeroRow = true;
            for (const auto& [j, v] : a.row(i)) {
                (void)j;
                if (std::abs(v) > eps) {
                    zeroRow = false;
                    break;
                }
            }
            if (zeroRow && std::abs(b[i]) > eps) {
                return true;
            }
        }
        return false;
    }

    static void applyPartialColumnPivoting(SparseMatrix& a, std::vector<double>& b, double eps) {
        const std::size_t n = a.size();
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t pivot = i;
            double best = std::abs(a.get(i, i));
            for (std::size_t r = i + 1; r < n; ++r) {
                const double cand = std::abs(a.get(r, i));
                if (cand > best) {
                    best = cand;
                    pivot = r;
                }
            }
            if (best <= eps) {
                continue;
            }
            if (pivot != i) {
                a.swapRows(i, pivot);
                std::swap(b[i], b[pivot]);
            }
        }
    }

    static RankResult classifyByRanks(const SparseMatrix& a, const std::vector<double>& b, double eps) {
        const std::size_t n = a.size();
        std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1, 0.0));

        for (std::size_t i = 0; i < n; ++i) {
            for (const auto& [j, v] : a.row(i)) {
                aug[i][j] = v;
            }
            aug[i][n] = b[i];
        }

        std::size_t row = 0;
        for (std::size_t col = 0; col < n && row < n; ++col) {
            std::size_t pivot = row;
            double mx = std::abs(aug[row][col]);
            for (std::size_t r = row + 1; r < n; ++r) {
                const double cur = std::abs(aug[r][col]);
                if (cur > mx) {
                    mx = cur;
                    pivot = r;
                }
            }
            if (mx <= eps) {
                continue;
            }
            if (pivot != row) {
                std::swap(aug[pivot], aug[row]);
            }
            for (std::size_t r = row + 1; r < n; ++r) {
                if (std::abs(aug[r][col]) <= eps) {
                    continue;
                }
                const double k = aug[r][col] / aug[row][col];
                for (std::size_t c = col; c <= n; ++c) {
                    aug[r][c] -= k * aug[row][c];
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

        RankResult out;
        if (rankA < rankAug) {
            out.type = SolutionType::NONE;
        } else if (rankA == rankAug && rankA < n) {
            out.type = SolutionType::INFINITE;
        } else if (rankA == n) {
            out.type = SolutionType::UNIQUE;
        }
        out.rankA = static_cast<int>(rankA);
        out.rankAb = static_cast<int>(rankAug);
        return out;
    }

    static SolveResult solveByJacobi(SparseMatrix a, std::vector<double> b, const JacobiConfig& cfg) {
        SolveResult res;
        const std::size_t n = a.size();

        if (n == 0 || b.size() != n) {
            res.message = "Некорректные размеры A и b.";
            return res;
        }

        if (cfg.usePivot) {
            applyPartialColumnPivoting(a, b, cfg.zeroEps);
        }

        if (hasInconsistentZeroRow(a, b, cfg.zeroEps)) {
            res.type = SolutionType::NONE;
            res.message = "Система несовместна (нулевая строка при ненулевом b).";
            return res;
        }

        if (isStrictlyDiagonallyDominant(a, cfg.zeroEps)) {
            res.type = SolutionType::UNIQUE;
            res.rankA = static_cast<int>(n);
            res.rankAb = static_cast<int>(n);
        } else if (n <= cfg.rankCheckLimit) {
            const RankResult rankRes = classifyByRanks(a, b, cfg.zeroEps);
            res.type = rankRes.type;
            res.rankA = rankRes.rankA;
            res.rankAb = rankRes.rankAb;
            if (res.type == SolutionType::NONE) {
                res.message = "rank(A) < rank([A|b]): решений нет.";
                return res;
            }
            if (res.type == SolutionType::INFINITE) {
                res.message = "rank(A)=rank([A|b])<n: бесконечно много решений.";
                return res;
            }
        }

        std::vector<double> prev(n, 0.0), next(n, 0.0);
        res.memoryBytes = a.memoryBytes() + b.capacity() * sizeof(double) + prev.capacity() * sizeof(double) +
                        next.capacity() * sizeof(double);

        double oldDelta = std::numeric_limits<double>::infinity();
        for (std::size_t it = 1; it <= cfg.maxIterations; ++it) {
            double delta = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                double diag = 0.0;
                double sum = 0.0;
                for (const auto& [j, v] : a.row(i)) {
                    if (j == i) {
                        diag = v;
                    } else {
                        sum += v * prev[j];
                    }
                }
                if (std::abs(diag) <= cfg.zeroEps) {
                    res.message = "На диагонали ноль/почти ноль, метод Якоби неприменим.";
                    res.iterations = it - 1;
                    res.achievedDelta = delta;
                    return res;
                }
                const double v = (b[i] - sum) / diag;
                if (!std::isfinite(v)) {
                    res.message = "Численная ошибка (NaN/Inf).";
                    res.iterations = it - 1;
                    return res;
                }
                next[i] = v;
                delta = std::max(delta, std::abs(next[i] - prev[i]));
            }

            res.iterations = it;
            res.achievedDelta = delta;
            if (delta < cfg.tolerance) {
                res.success = true;
                res.converged = true;
                res.x = next;
                res.residualInf = residualInfNorm(a, res.x, b);
                res.message = "Сходимость достигнута.";
                return res;
            }

            if (it > 8 && delta > oldDelta * 1e6) {
                res.message = "Наблюдается расходимость итераций.";
                res.x = next;
                res.residualInf = residualInfNorm(a, res.x, b);
                return res;
            }

            oldDelta = delta;
            std::swap(prev, next);
        }

        res.x = prev;
        res.residualInf = residualInfNorm(a, res.x, b);
        res.message = "Достигнут предел итераций без сходимости.";
        return res;
    }

    static SparseMatrix buildTridiagonal(std::size_t n, double diag, double off) {
        SparseMatrix a(n);
        for (std::size_t i = 0; i < n; ++i) {
            a.set(i, i, diag);
            if (i > 0) {
                a.set(i, i - 1, off);
            }
            if (i + 1 < n) {
                a.set(i, i + 1, off);
            }
        }
        return a;
    }

    static SparseMatrix buildHilbert(std::size_t n) {
        SparseMatrix h(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                h.set(i, j, 1.0 / static_cast<double>(i + j + 1));
            }
        }
        return h;
    }

    static std::vector<double> randomVector(std::size_t n, std::mt19937_64& rng) {
        std::uniform_real_distribution<double> dist(-10.0, 10.0);
        std::vector<double> x(n, 0.0);
        for (double& v : x) {
            v = dist(rng);
        }
        return x;
    }

    static SparseMatrix randomBandedDominantMatrix(std::size_t n, std::mt19937_64& rng) {
        std::uniform_real_distribution<double> dist(-2.0, 2.0);
        SparseMatrix a(n);
        for (std::size_t i = 0; i < n; ++i) {
            double offAbs = 0.0;
            if (i > 0) {
                const double v = dist(rng);
                a.set(i, i - 1, v);
                offAbs += std::abs(v);
            }
            if (i + 1 < n) {
                const double v = dist(rng);
                a.set(i, i + 1, v);
                offAbs += std::abs(v);
            }
            if (i > 1) {
                const double v = dist(rng) * 0.4;
                a.set(i, i - 2, v);
                offAbs += std::abs(v);
            }
            if (i + 2 < n) {
                const double v = dist(rng) * 0.4;
                a.set(i, i + 2, v);
                offAbs += std::abs(v);
            }
            a.set(i, i, offAbs + 5.0);
        }
        return a;
    }

    static std::string statusToString(SolutionType t) {
        switch (t) {
            case SolutionType::UNIQUE:
                return "Единственное решение";
            case SolutionType::INFINITE:
                return "Бесконечно много решений";
            case SolutionType::NONE:
                return "Решений нет";
            default:
                return "Не удалось строго классифицировать";
        }
    }

    static void printRankInfo(const SolveResult& result) {
        if (result.rankA >= 0 && result.rankAb >= 0) {
            std::cout << "rank(A) = " << result.rankA << ", rank([A|b]) = " << result.rankAb << "\n";
        } else {
            std::cout << "rank(A) / rank([A|b]) не вычислялись (большой размер или быстрый режим).\n";
        }
    }

    static void runSingleRandomCase() {
        std::size_t n;
        std::cout << "Введите размер матрицы N: ";
        std::cin >> n;
        if (n < 2) {
            std::cout << "N должно быть >= 2\n";
            return;
        }

        std::mt19937_64 rng(std::random_device{}());
        SparseMatrix a = randomBandedDominantMatrix(n, rng);
        std::vector<double> xRef = randomVector(n, rng);
        std::vector<double> b = multiply(a, xRef);

        JacobiConfig cfg;
        cfg.tolerance = 1e-10;
        cfg.maxIterations = 100000;
        cfg.rankCheckLimit = 200;

        const auto start = std::chrono::high_resolution_clock::now();
        SolveResult result = solveByJacobi(a, b, cfg);
        const auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Статус: " << statusToString(result.type) << "\n";
        printRankInfo(result);
        std::cout << "Время: " << std::chrono::duration<double>(end - start).count() << " сек\n";
        if (result.success) {
            std::cout << "e = ||x_ref - x_star||_2 = " << formatScientific(l2ErrorNorm(xRef, result.x)) << "\n";
            std::cout << "Невязка ||A*x-b||inf = " << result.residualInf << "\n";
        } else {
            std::cout << "Сообщение: " << result.message << "\n";
        }
    }

    static void runBasicTests() {
        std::cout << "\n=== БАЗОВЫЕ ТЕСТЫ ===\n";

        {
            SparseMatrix a(3);
            a.set(0, 0, 10); a.set(0, 1, 1); a.set(0, 2, 1);
            a.set(1, 0, 2);  a.set(1, 1, 10); a.set(1, 2, 1);
            a.set(2, 0, 2);  a.set(2, 1, 2);  a.set(2, 2, 10);
            std::vector<double> xRef = {1.0, 1.0, 1.0};
            std::vector<double> b = multiply(a, xRef);

            SolveResult r = solveByJacobi(a, b, JacobiConfig{});
            std::cout << "[Тест 1] Классическая 3x3 система: " << statusToString(r.type);
            if (r.success) {
                std::cout << ", e = " << formatScientific(l2ErrorNorm(xRef, r.x));
            }
            std::cout << "\n";
        }

        {
            SparseMatrix a(2);
            a.set(0, 0, 1); a.set(0, 1, 1);
            a.set(1, 0, 2); a.set(1, 1, 2);
            std::vector<double> b = {2.0, 4.0};

            SolveResult r = solveByJacobi(a, b, JacobiConfig{});
            std::cout << "[Тест 2] Вырожденная совместная: " << statusToString(r.type) << "\n";
        }

        {
            SparseMatrix a(2);
            a.set(0, 0, 1); a.set(0, 1, 1);
            a.set(1, 0, 2); a.set(1, 1, 2);
            std::vector<double> b = {2.0, 5.0};

            SolveResult r = solveByJacobi(a, b, JacobiConfig{});
            std::cout << "[Тест 3] Вырожденная несовместная: " << statusToString(r.type) << "\n";
        }

        {
            const std::size_t n = 8;
            std::mt19937_64 rng(42);
            SparseMatrix a = randomBandedDominantMatrix(n, rng);
            std::vector<double> xRef = randomVector(n, rng);
            std::vector<double> b = multiply(a, xRef);

            SolveResult r = solveByJacobi(a, b, JacobiConfig{});
            std::cout << "[Тест 4] Случайная невырожденная N=8: " << statusToString(r.type);
            if (r.success) {
                std::cout << ", e = " << formatScientific(l2ErrorNorm(xRef, r.x));
            }
            std::cout << "\n";
        }
    }

    static void runStressTests() {
        const std::vector<std::size_t> sizes = {1000, 10000};

        for (const std::size_t n : sizes) {
            std::cout << "\nN = " << n << "\n";
            try {
                const auto createStart = std::chrono::high_resolution_clock::now();
                SparseMatrix a = buildTridiagonal(n, 10.0, -1.0);
                std::vector<double> xRef(n, 1.0);
                std::vector<double> b = multiply(a, xRef);
                const auto createEnd = std::chrono::high_resolution_clock::now();

                JacobiConfig cfg;
                cfg.usePivot = false;
                cfg.tolerance = (n == 1000) ? 1e-8 : 1e-6;
                cfg.maxIterations = (n == 1000) ? 5000 : 3000;
                cfg.rankCheckLimit = 80;

                const auto solveStart = std::chrono::high_resolution_clock::now();
                SolveResult result = solveByJacobi(std::move(a), std::move(b), cfg);
                const auto solveEnd = std::chrono::high_resolution_clock::now();

                std::cout << "Статус: " << statusToString(result.type) << "\n";
                std::cout << "Подготовка данных: "
                        << std::chrono::duration<double>(createEnd - createStart).count() << " сек\n";
                std::cout << "Решение: " << std::chrono::duration<double>(solveEnd - solveStart).count()
                        << " сек\n";
                std::cout << "Память (оценка): " << result.memoryBytes << " байт\n";

                if (result.success) {
                    std::cout << "e = ||x_ref - x_star||_2 = " << formatScientific(l2ErrorNorm(xRef, result.x)) << "\n";
                    std::cout << "Невязка ||A*x-b||inf = " << result.residualInf << "\n";
                } else {
                    std::cout << "Сообщение: " << result.message << "\n";
                }
            } catch (const std::bad_alloc&) {
                std::cout << "Недостаточно памяти для матрицы " << n << "x" << n << "\n";
            } catch (const std::exception& e) {
                std::cout << "Ошибка: " << e.what() << "\n";
            }
        }
    }

    static void runHilbertTest() {
        std::size_t n;
        std::cout << "Введите размер матрицы Гильберта N: ";
        std::cin >> n;
        if (n < 2) {
            std::cout << "N должно быть >= 2\n";
            return;
        }

        std::mt19937_64 rng(123456);
        SparseMatrix a = buildHilbert(n);
        std::vector<double> xRef = randomVector(n, rng);
        std::vector<double> b = multiply(a, xRef);

        JacobiConfig cfg;
        cfg.tolerance = 1e-12;
        cfg.maxIterations = 5000;
        cfg.rankCheckLimit = 200;

        const auto start = std::chrono::high_resolution_clock::now();
        SolveResult result = solveByJacobi(a, b, cfg);
        const auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Статус: " << statusToString(result.type) << "\n";
        std::cout << "Время: " << std::chrono::duration<double>(end - start).count() << " сек\n";
        if (result.success) {
            std::cout << "e = ||x_ref - x_star||_2 = " << formatScientific(l2ErrorNorm(xRef, result.x)) << "\n";
        } else {
            std::cout << "Сообщение: " << result.message << "\n";
            std::cout << "Для матрицы Гильберта метод Якоби часто сходится медленно или расходится.\n";
        }
    }

    int main() {
        std::setlocale(LC_ALL, "Russian");
        std::cout << std::fixed << std::setprecision(10);

        std::cout << "Лабораторная работа: решение СЛАУ методом Якоби (частичный выбор)\n";
        std::cout << "1 - Один случайный запуск\n";
        std::cout << "2 - Базовые тесты\n";
        std::cout << "3 - Стресс-тесты (N=1000 и N=10000)\n";
        std::cout << "4 - Тест на матрице Гильберта\n";
        std::cout << "Ваш выбор: ";

        int choice = 0;
        std::cin >> choice;

        switch (choice) {
            case 1:
                runSingleRandomCase();
                break;
            case 2:
                runBasicTests();
                break;
            case 3:
                runStressTests();
                break;
            case 4:
                runHilbertTest();
                break;
            default:
                std::cout << "Неверный пункт меню\n";
                break;
        }

        return 0;
    }
