// Numerical integration (composite rectangles, trapezoids, Simpson)
// Variant 19: integral from 0 to pi of x^4 * exp(-x^2) dx

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using ld = long double;

constexpr ld PI = acosl(-1.0L);

ld f(ld x) {
	ld x2 = x * x;
	return x2 * x2 * expl(-x2);
}

ld rectangleFormula(ld a, ld b, int n) {
	ld h = (b - a) / static_cast<ld>(n);
	ld sum = 0.0L;
	for (int i = 0; i < n; ++i) {
		ld x = a + (static_cast<ld>(i) + 0.5L) * h;
		sum += f(x);
	}
	return h * sum;
}

ld trapezoidFormula(ld a, ld b, int n) {
	ld h = (b - a) / static_cast<ld>(n);
	ld sum = (f(a) + f(b)) / 2.0L;
	for (int i = 1; i < n; ++i) {
		ld x = a + static_cast<ld>(i) * h;
		sum += f(x);
	}
	return h * sum;
}

ld simpsonFormula(ld a, ld b, int n) {
	if (n % 2 != 0) {
		++n;
	}
	ld h = (b - a) / static_cast<ld>(n);
	ld sum = f(a) + f(b);
	for (int i = 1; i < n; ++i) {
		ld x = a + static_cast<ld>(i) * h;
		sum += (i % 2 == 0) ? 2.0L * f(x) : 4.0L * f(x);
	}
	return (h / 3.0L) * sum;
}

struct Result {
	std::string methodName;
	ld eps;
	int n;
	ld h;
	ld integral;
	ld rungeIntegral;
	ld error;
	int iterations;
};

Result rungeRule(
	const std::string& methodName,
	const std::function<ld(ld, ld, int)>& formula,
	int order,
	ld a,
	ld b,
	ld eps
) {
	int n = 2;
	ld oldIntegral = formula(a, b, n);
	int iterations = 0;
	const int maxIterations = 60;

	while (iterations < maxIterations) {
		n *= 2;
		++iterations;

		ld newIntegral = formula(a, b, n);
		ld denominator = std::pow(2.0L, static_cast<ld>(order)) - 1.0L;
		ld error = fabsl(newIntegral - oldIntegral) / denominator;
		ld rungeIntegral = newIntegral + (newIntegral - oldIntegral) / denominator;

		if (error <= eps) {
			return Result{
				methodName,
				eps,
				n,
				(b - a) / static_cast<ld>(n),
				newIntegral,
				rungeIntegral,
				error,
				iterations
			};
		}

		oldIntegral = newIntegral;
	}

	return Result{
		methodName,
		eps,
		n,
		(b - a) / static_cast<ld>(n),
		oldIntegral,
		oldIntegral,
		std::numeric_limits<ld>::quiet_NaN(),
		iterations
	};
}

void printResult(const Result& result) {
	std::cout << "Метод: " << result.methodName << '\n';
	std::cout << "eps = " << result.eps << '\n';
	std::cout << "n = " << result.n << '\n';
	std::cout << "h = " << result.h << '\n';
	std::cout << "I по формуле = " << result.integral << '\n';
	std::cout << "I с уточнением по Рунге = " << result.rungeIntegral << '\n';
	std::cout << "Оценка погрешности по Рунге = " << result.error << '\n';
	std::cout << "Количество удвоений n = " << result.iterations << '\n';
	std::cout << "------------------------------------------------------------\n";
}

int main() {
	ld a = 0.0L;
	ld b = PI;

	std::vector<ld> epsValues = {0.01L, 0.001L, 0.0001L, 0.00001L, 0.000001L, 0.000001L, 0.00000001L};

	std::cout << std::fixed << std::setprecision(12);
	std::cout << "Лабораторная работа №7\n";
	std::cout << "Численное интегрирование (вариант 19)\n";
	std::cout << "Интеграл: x^4 * exp(-x^2) на [0, pi]\n";
	std::cout << "------------------------------------------------------------\n\n";

	for (ld eps : epsValues) {
		std::cout << "Требуемая точность eps = " << eps << '\n';
		std::cout << "============================================================\n";

		Result rectangleResult = rungeRule(
			"Составная формула прямоугольников (середины)",
			rectangleFormula,
			2,
			a,
			b,
			eps
		);

		Result trapezoidResult = rungeRule(
			"Составная формула трапеций",
			trapezoidFormula,
			2,
			a,
			b,
			eps
		);

		Result simpsonResult = rungeRule(
			"Составная формула Симпсона",
			simpsonFormula,
			4,
			a,
			b,
			eps
		);

		printResult(rectangleResult);
		printResult(trapezoidResult);
		printResult(simpsonResult);

		std::cout << '\n';
	}

	return 0;
}
