import argparse
import numpy as np
import scipy.linalg


def parse_args():
    parser = argparse.ArgumentParser(
        description="Решение СЛАУ методом Гаусса с частичным выбором главного элемента"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    research_parser = subparsers.add_parser("research", help="Провести серию экспериментов")
    research_parser.add_argument("--n", type=int, required=True, help="Размер матрицы")
    research_parser.add_argument(
        "--type",
        choices=["random", "hilbert", "diagonal-dominant"],
        default="random",
        help="Тип матрицы"
    )
    research_parser.add_argument("--interval", type=float, default=100.0, help="Верхняя граница генерации")
    research_parser.add_argument("--seed", type=int, default=1, help="Начальное значение генератора")

    solve_parser = subparsers.add_parser("solve", help="Решить одну систему")
    solve_parser.add_argument("--n", type=int, required=True, help="Размер матрицы")
    solve_parser.add_argument(
        "--type",
        choices=["random", "hilbert", "diagonal-dominant"],
        default="random",
        help="Тип матрицы"
    )
    solve_parser.add_argument("--interval", type=float, default=100.0, help="Верхняя граница генерации")
    solve_parser.add_argument("--seed", type=int, default=1, help="Начальное значение генератора")

    return parser.parse_args()


def validate_args(args):
    if args.n < 2:
        raise ValueError("Размер матрицы должен быть не меньше 2")
    if args.interval <= 0:
        raise ValueError("Интервал генерации должен быть положительным")


def generate_random_system(n, interval=100.0, seed=1):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0, interval, size=(n, n))
    b = rng.uniform(0, interval, size=n)
    return a, b


def generate_hilbert_system(n):
    idx = np.arange(1, n + 1)
    a = 1.0 / (idx[:, None] + idx[None, :] - 1)
    x_exact = np.ones(n)
    b = a @ x_exact
    return a, b, x_exact


def generate_diagonal_dominant_system(n, interval=100.0, seed=1):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0, interval, size=(n, n))
    for i in range(n):
        row_sum = np.sum(np.abs(a[i])) - abs(a[i, i])
        a[i, i] = row_sum + rng.uniform(1.0, interval)
    b = rng.uniform(0, interval, size=n)
    return a, b


def gauss_partial_pivot(a, b):
    a = a.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)

    for k in range(n):
        pivot_row = k + np.argmax(np.abs(a[k:, k]))
        if abs(a[pivot_row, k]) < 1e-15:
            raise ValueError("Матрица вырождена или близка к вырожденной")

        if pivot_row != k:
            a[[k, pivot_row]] = a[[pivot_row, k]]
            b[[k, pivot_row]] = b[[pivot_row, k]]

        pivot = a[k, k]
        a[k, k:] /= pivot
        b[k] /= pivot

        if k + 1 < n:
            factors = a[k + 1:, k].copy()
            a[k + 1:, k:] -= factors[:, None] * a[k, k:]
            b[k + 1:] -= factors * b[k]

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = b[i] - a[i, i + 1:] @ x[i + 1:]
    return x


def relative_residual(a, x, b):
    numerator = np.linalg.norm(a @ x - b, ord=2)
    denominator = np.linalg.norm(b, ord=2)
    return numerator / denominator if denominator != 0 else numerator


def relative_difference(x1, x2):
    return np.linalg.norm(x1 - x2, ord=2) / np.linalg.norm(x2, ord=2)


def analyze_system(a, b):
    x_gauss = gauss_partial_pivot(a, b)
    x_scipy = scipy.linalg.solve(a, b)
    return {
        "gauss_solution": x_gauss,
        "scipy_solution": x_scipy,
        "residual": relative_residual(a, x_gauss, b),
        "condition_number": np.linalg.cond(a),
        "difference": relative_difference(x_gauss, x_scipy),
    }


def make_system(matrix_type, n, interval, seed):
    if matrix_type == "random":
        return generate_random_system(n, interval, seed), None
    if matrix_type == "diagonal-dominant":
        return generate_diagonal_dominant_system(n, interval, seed), None
    if matrix_type == "hilbert":
        a, b, x_exact = generate_hilbert_system(n)
        return (a, b), x_exact
    raise ValueError("Неизвестный тип матрицы")


def run_solve(args):
    (a, b), x_exact = make_system(args.type, args.n, args.interval, args.seed)
    result = analyze_system(a, b)

    print(f"Тип матрицы: {args.type}")
    print(f"Размер: {args.n}")
    print("-" * 70)
    print(f"{'i':<5}{'Метод Гаусса':<24}{'SciPy':<24}")
    print("-" * 70)
    for i, (g, s) in enumerate(zip(result["gauss_solution"], result["scipy_solution"]), start=1):
        print(f"{i:<5}{g:<24.10f}{s:<24.10f}")
    print("-" * 70)
    print(f"Относительная невязка: {result['residual']:.6e}")
    print(f"Число обусловленности: {result['condition_number']:.6e}")
    print(f"Относительное отличие от SciPy: {result['difference']:.6e}")

    if x_exact is not None:
        max_error = np.max(np.abs(result["gauss_solution"] - x_exact))
        print(f"Максимальная ошибка относительно точного решения: {max_error:.6e}")


def run_research(args):
    sizes = [args.n]
    if args.type == "hilbert":
        print(f"{'n':<8}{'Невязка':<18}{'cond(A)':<18}{'Отличие от SciPy':<20}{'Макс. ошибка':<18}")
        print("-" * 82)
        for n in sizes:
            (a, b), x_exact = make_system(args.type, n, args.interval, args.seed)
            result = analyze_system(a, b)
            max_error = np.max(np.abs(result["gauss_solution"] - x_exact))
            print(
                f"{n:<8}{result['residual']:<18.6e}{result['condition_number']:<18.6e}"
                f"{result['difference']:<20.6e}{max_error:<18.6e}"
            )
    else:
        print(f"{'n':<8}{'Невязка':<18}{'cond(A)':<18}{'Отличие от SciPy':<20}")
        print("-" * 64)
        for n in sizes:
            (a, b), _ = make_system(args.type, n, args.interval, args.seed)
            result = analyze_system(a, b)
            print(
                f"{n:<8}{result['residual']:<18.6e}{result['condition_number']:<18.6e}"
                f"{result['difference']:<20.6e}"
            )


def main():
    args = parse_args()
    validate_args(args)
    if args.command == "solve":
        run_solve(args)
    else:
        run_research(args)


if __name__ == "__main__":
    main()