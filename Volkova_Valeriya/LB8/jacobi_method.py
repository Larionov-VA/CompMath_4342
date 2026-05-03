import argparse
import numpy as np
import scipy.linalg


def parse_args():
    parser = argparse.ArgumentParser(
        description="Решение СЛАУ методом Якоби с исследованием сходимости"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    research_parser = subparsers.add_parser("research", help="Провести серию экспериментов")
    research_parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=[10, 100, 300],
        help="Список размеров матриц для исследования"
    )
    research_parser.add_argument(
        "--type",
        choices=["all", "random", "hilbert", "diagonal-dominant"],
        default="all",
        help="Тип матрицы"
    )
    research_parser.add_argument("--interval", type=float, default=100.0, help="Верхняя граница генерации")
    research_parser.add_argument("--seed", type=int, default=1, help="Начальное значение генератора")
    research_parser.add_argument("--tol", type=float, default=1e-10, help="Точность остановки")
    research_parser.add_argument("--max-iter", type=int, default=20000, help="Максимум итераций")

    solve_parser = subparsers.add_parser("solve", help="Решить одну систему")
    solve_parser.add_argument("--n", type=int, required=True, help="Размер матрицы")
    solve_parser.add_argument(
        "--type",
        choices=["random", "hilbert", "diagonal-dominant"],
        default="diagonal-dominant",
        help="Тип матрицы"
    )
    solve_parser.add_argument("--interval", type=float, default=100.0, help="Верхняя граница генерации")
    solve_parser.add_argument("--seed", type=int, default=1, help="Начальное значение генератора")
    solve_parser.add_argument("--tol", type=float, default=1e-10, help="Точность остановки")
    solve_parser.add_argument("--max-iter", type=int, default=20000, help="Максимум итераций")

    return parser.parse_args()


def validate_args(args):
    if getattr(args, "n", 2) < 2:
        raise ValueError("Размер матрицы должен быть не меньше 2")
    sizes = getattr(args, "sizes", None)
    if sizes is not None and any(size < 2 for size in sizes):
        raise ValueError("Все размеры матриц должны быть не меньше 2")
    if args.interval <= 0:
        raise ValueError("Интервал генерации должен быть положительным")
    if args.tol <= 0:
        raise ValueError("Точность должна быть положительной")
    if args.max_iter <= 0:
        raise ValueError("Максимум итераций должен быть положительным")


def generate_random_system(n, interval=100.0, seed=1):
    rng = np.random.default_rng(seed)
    a = rng.uniform(-interval, interval, size=(n, n))
    b = rng.uniform(-interval, interval, size=n)
    return a, b


def generate_hilbert_system(n):
    idx = np.arange(1, n + 1)
    a = 1.0 / (idx[:, None] + idx[None, :] - 1)
    x_exact = np.ones(n)
    b = a @ x_exact
    return a, b, x_exact


def generate_diagonal_dominant_system(n, interval=100.0, seed=1):
    rng = np.random.default_rng(seed)
    a = rng.uniform(-interval, interval, size=(n, n))

    for i in range(n):
        row_sum = np.sum(np.abs(a[i])) - abs(a[i, i])
        a[i, i] = row_sum + rng.uniform(1.0, interval)

    b = rng.uniform(-interval, interval, size=n)
    return a, b


def relative_residual(a, x, b):
    numerator = np.linalg.norm(a @ x - b, ord=2)
    denominator = np.linalg.norm(b, ord=2)
    return numerator / denominator if denominator != 0 else numerator


def relative_difference(x1, x2):
    denominator = np.linalg.norm(x2, ord=2)
    if denominator == 0:
        return np.linalg.norm(x1 - x2, ord=2)
    return np.linalg.norm(x1 - x2, ord=2) / denominator


def iteration_matrix_inf_norm(a):
    d = np.diag(a)
    if np.any(np.abs(d) < 1e-15):
        return float("inf")
    off_diag = a - np.diag(d)
    b_matrix = -((off_diag.T) / d).T
    return np.linalg.norm(b_matrix, ord=np.inf)


def jacobi_method(a, b, tol=1e-10, max_iter=20000):
    a = a.astype(float).copy()
    b = b.astype(float).copy()
    d = np.diag(a)

    if np.any(np.abs(d) < 1e-15):
        raise ValueError("На диагонали есть нулевой элемент, метод Якоби неприменим")

    r = a - np.diagflat(d)
    x = np.zeros_like(b)
    previous_diff = np.inf

    for iteration in range(1, max_iter + 1):
        x_new = (b - r @ x) / d

        if not np.all(np.isfinite(x_new)):
            return x, iteration, False, float("inf")

        diff = np.linalg.norm(x_new - x, ord=np.inf)
        if diff < tol:
            return x_new, iteration, True, diff

        if iteration > 25 and (diff > previous_diff * 10 or np.linalg.norm(x_new, ord=np.inf) > 1e100):
            return x_new, iteration, False, diff

        previous_diff = diff
        x = x_new

    return x, max_iter, False, diff


def make_system(matrix_type, n, interval, seed):
    if matrix_type == "random":
        return generate_random_system(n, interval, seed), None
    if matrix_type == "diagonal-dominant":
        return generate_diagonal_dominant_system(n, interval, seed), None
    if matrix_type == "hilbert":
        a, b, x_exact = generate_hilbert_system(n)
        return (a, b), x_exact
    raise ValueError("Неизвестный тип матрицы")


def analyze_system(a, b, tol=1e-10, max_iter=20000):
    q_inf = iteration_matrix_inf_norm(a)
    x_jacobi, iterations, converged, final_diff = jacobi_method(a, b, tol, max_iter)
    condition_number = np.linalg.cond(a)

    result = {
        "jacobi_solution": x_jacobi,
        "iterations": iterations,
        "converged": converged,
        "final_diff": final_diff,
        "iteration_matrix_norm": q_inf,
        "condition_number": condition_number,
    }

    if converged:
        x_scipy = scipy.linalg.solve(a, b)
        result["scipy_solution"] = x_scipy
        result["residual"] = relative_residual(a, x_jacobi, b)
        result["difference"] = relative_difference(x_jacobi, x_scipy)
    else:
        result["scipy_solution"] = None
        result["residual"] = None
        result["difference"] = None

    return result


def print_solution_table(result, n, x_exact=None):
    limit = min(n, 20)
    print(f"{'i':<5}{'Метод Якоби':<24}{'SciPy':<24}")
    print("-" * 70)
    if result["converged"]:
        for i in range(limit):
            print(f"{i + 1:<5}{result['jacobi_solution'][i]:<24.10f}{result['scipy_solution'][i]:<24.10f}")
        if n > limit:
            print("... показаны только первые 20 компонент ...")
    else:
        print("Метод Якоби не сошелся за заданное число итераций.")
    print("-" * 70)
    print(f"Сходимость: {'да' if result['converged'] else 'нет'}")
    print(f"Число итераций: {result['iterations']}")
    print(f"||B||_inf: {result['iteration_matrix_norm']:.6e}")
    print(f"Число обусловленности: {result['condition_number']:.6e}")
    if result["residual"] is not None:
        print(f"Относительная невязка: {result['residual']:.6e}")
        print(f"Относительное отличие от SciPy: {result['difference']:.6e}")
    if x_exact is not None and result["converged"]:
        max_error = np.max(np.abs(result["jacobi_solution"] - x_exact))
        print(f"Максимальная ошибка относительно точного решения: {max_error:.6e}")


def run_solve(args):
    (a, b), x_exact = make_system(args.type, args.n, args.interval, args.seed)
    result = analyze_system(a, b, args.tol, args.max_iter)

    print(f"Тип матрицы: {args.type}")
    print(f"Размер: {args.n}")
    print_solution_table(result, args.n, x_exact)


def run_research(args):
    matrix_types = [args.type] if args.type != "all" else ["random", "hilbert", "diagonal-dominant"]

    print(
        f"{'n':<8}{'Тип':<22}{'Сошелся':<10}{'Итерации':<12}{'||B||_inf':<16}"
        f"{'Невязка':<16}{'cond(A)':<16}{'Отличие':<16}"
    )
    print("-" * 116)

    for n in args.sizes:
        for matrix_type in matrix_types:
            (a, b), _ = make_system(matrix_type, n, args.interval, args.seed)
            result = analyze_system(a, b, args.tol, args.max_iter)
            residual_text = f"{result['residual']:.3e}" if result["residual"] is not None else "-"
            difference_text = f"{result['difference']:.3e}" if result["difference"] is not None else "-"
            print(
                f"{n:<8}{matrix_type:<22}{('да' if result['converged'] else 'нет'):<10}"
                f"{result['iterations']:<12}{result['iteration_matrix_norm']:<16.3e}"
                f"{residual_text:<16}{result['condition_number']:<16.3e}{difference_text:<16}"
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
