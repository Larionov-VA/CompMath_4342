import math
import sys
from dataclasses import dataclass
from typing import Callable, List

"""

    Тема:
    Численное интегрирование.
    Составные квадратурные формулы прямоугольников, трапеций и Симпсона.
    Оценка погрешности по правилу Рунге.

    Интеграл от 0 до pi:
        x^2 * exp(-x^2) dx

    Точность задается через запуск программы.

    Примеры запуска:
        python main.py 0.01
        python main.py 0.01 0.001 0.0001
"""


# Константа pi.
PI = math.pi


def f(x: float) -> float:
    """
    Подынтегральная функция.

    По условию:
        f(x) = x^2 * exp(-x^2)

    На отрезке [0; pi] функция определена и непрерывна,
    поэтому для нее можно применять формулы прямоугольников,
    трапеций и Симпсона.
    """

    return x ** 2 * math.exp(-x ** 2)


def rectangle_formula(a: float, b: float, n: int) -> float:
    """
    Составная формула прямоугольников.

    Используется формула средних прямоугольников:
        I ≈ h * сумма f(a + (i + 0.5) * h), i = 0, 1, ..., n - 1

    Здесь:
        a, b — границы интегрирования;
        n — количество частей разбиения;
        h = (b - a) / n — шаг сетки.
    """

    h = (b - a) / n
    total = 0.0

    for i in range(n):
        x = a + (i + 0.5) * h
        total += f(x)

    return h * total


def trapezoid_formula(a: float, b: float, n: int) -> float:
    """
    Составная формула трапеций.

    Формула:
        I ≈ h * ((f(a) + f(b)) / 2 + сумма f(a + i*h), i = 1, 2, ..., n - 1)

    Идея:
    На каждом маленьком отрезке график функции заменяется прямой линией,
    поэтому площадь считается как площадь трапеции.
    """

    h = (b - a) / n
    total = (f(a) + f(b)) / 2.0

    for i in range(1, n):
        x = a + i * h
        total += f(x)

    return h * total


def simpson_formula(a: float, b: float, n: int) -> float:
    """
    Составная формула Симпсона.

    Формула:
        I ≈ h/3 * (f(x0) + 4*f(x1) + 2*f(x2) + ... + 4*f(x(n-1)) + f(xn))

    Важно:
    Для формулы Симпсона количество отрезков n обязательно должно быть четным.
    """

    if n % 2 != 0:
        raise ValueError("Для формулы Симпсона n должно быть четным")

    h = (b - a) / n
    total = f(a) + f(b)

    for i in range(1, n):
        x = a + i * h

        if i % 2 == 1:
            total += 4.0 * f(x)
        else:
            total += 2.0 * f(x)

    return h / 3.0 * total


@dataclass
class Method:
    """
    Описание одного метода интегрирования.

    Поля:
        name — название метода;
        formula — функция для вычисления интеграла;
        order — порядок точности метода;
        runge_sign — знак уточнения по Рунге.

    По заданию:
        для прямоугольников итоговое значение: I_h + погрешность;
        для трапеций итоговое значение: I_h - погрешность;
        для Симпсона итоговое значение: I_h - погрешность.
    """

    name: str
    formula: Callable[[float, float, int], float]
    order: int
    runge_sign: int


@dataclass
class Result:
    """
    Результат работы одного метода.

    Поля:
        method_name — название метода;
        eps — заданная точность;
        n_coarse — число отрезков на грубой сетке;
        n_fine — число отрезков на мелкой сетке;
        h — шаг мелкой сетки;
        i_2h — значение интеграла на грубой сетке;
        i_h — значение интеграла на мелкой сетке;
        error — оценка погрешности по правилу Рунге;
        runge_integral — итоговое уточненное значение интеграла;
        checks — сколько раз проверяли n.
    """

    method_name: str
    eps: float
    n_coarse: int
    n_fine: int
    h: float
    i_2h: float
    i_h: float
    error: float
    runge_integral: float
    checks: int


def parse_eps_values(args: List[str]) -> List[float]:
    """
    Получает точности из командной строки.

    Пример:
        python main.py 0.01 0.001 0.0001
    """

    if not args:
        print("Ошибка: точность нужно задать через запуск программы.")
        print()
        print("Примеры запуска:")
        print("  python main.py 0.01")
        print("  python main.py 0.01 0.001 0.0001")
        sys.exit(1)

    eps_values = []

    for value in args:
        try:
            # replace нужен, чтобы можно было ввести и 0.001, и 0,001.
            eps = float(value.replace(",", "."))
        except ValueError:
            print(f"Ошибка: значение '{value}' не является числом.")
            sys.exit(1)

        if eps <= 0:
            print(f"Ошибка: точность должна быть положительной, получено {value}.")
            sys.exit(1)

        eps_values.append(eps)

    return eps_values


def runge_rule(method: Method, a: float, b: float, eps: float) -> Result:
    """
    Подбор n по правилу Рунге.

    Для каждого n считаются два значения:
        I_2h — значение на более грубой сетке;
        I_h  — значение на более мелкой сетке.

    Оценка погрешности:
        error = |I_h - I_2h| / (2^p - 1)

    где p — порядок точности метода:
        p = 2 для прямоугольников и трапеций;
        p = 4 для Симпсона.

    Цикл продолжается до тех пор, пока:
        error <= eps
    """

    denominator = 2 ** method.order - 1

    # Начинаем с n = 4.
    # Тогда грубая сетка имеет n/2 = 2 отрезка,
    # а для Симпсона оба значения считаются при четном количестве отрезков.
    n_fine = 4
    checks = 0

    while True:
        n_coarse = n_fine // 2

        i_2h = method.formula(a, b, n_coarse)
        i_h = method.formula(a, b, n_fine)

        error = abs(i_h - i_2h) / denominator
        runge_integral = i_h + method.runge_sign * error

        checks += 1

        if error <= eps:
            return Result(
                method_name=method.name,
                eps=eps,
                n_coarse=n_coarse,
                n_fine=n_fine,
                h=(b - a) / n_fine,
                i_2h=i_2h,
                i_h=i_h,
                error=error,
                runge_integral=runge_integral,
                checks=checks
            )

        # Берем степени двойки
        n_fine *= 2


def print_header(eps_values: List[float]) -> None:
    """
    Печатает общую информацию о задаче.
    """

    print("=" * 118)
    print("Интеграл: ∫[0; pi] x^2 * exp(-x^2) dx")
    print("Точности из запуска:", ", ".join(f"{eps:g}" for eps in eps_values))
    print("=" * 118)
    print()
    print("Обозначения в таблице:")
    print("  n/2       — число отрезков на грубой сетке")
    print("  n         — число отрезков на мелкой сетке")
    print("  h         — шаг мелкой сетки")
    print("  I_2h      — значение интеграла на грубой сетке")
    print("  I_h       — значение интеграла на мелкой сетке")
    print("  Погр.     — оценка погрешности по правилу Рунге")
    print("  Итоговый I — уточненное значение интеграла")
    print()


def print_results_table(eps: float, results: List[Result]) -> None:
    """
    Печатает таблицу результатов для одной точности eps.
    """

    print("=" * 118)
    print(f"Требуемая точность eps = {eps:g}")
    print("-" * 118)

    print(
        f"{'Метод':<20}"
        f"{'n/2':>8}"
        f"{'n':>8}"
        f"{'h':>15}"
        f"{'I_2h':>18}"
        f"{'I_h':>18}"
        f"{'Погр.':>18}"
        f"{'Итоговый I':>18}"
    )

    print("-" * 118)

    for result in results:
        print(
            f"{result.method_name:<20}"
            f"{result.n_coarse:>8}"
            f"{result.n_fine:>8}"
            f"{result.h:>15.10f}"
            f"{result.i_2h:>18.10f}"
            f"{result.i_h:>18.10f}"
            f"{result.error:>18.10f}"
            f"{result.runge_integral:>18.10f}"
        )

    print("-" * 118)

    best = min(results, key=lambda item: item.n_fine)
    print(
        f"Меньше всего отрезков потребовалось методу: "
        f"{best.method_name}, n = {best.n_fine}."
    )
    print()


def main() -> None:
    """
    Главная функция программы.
    """

    eps_values = parse_eps_values(sys.argv[1:])

    # Границы интегрирования по варианту 18.
    a = 0.0
    b = PI

    methods = [
        Method("Прямоугольники", rectangle_formula, 2, 1),
        Method("Трапеции", trapezoid_formula, 2, -1),
        Method("Симпсон", simpson_formula, 4, -1),
    ]

    print_header(eps_values)

    for eps in eps_values:
        results = []

        for method in methods:
            result = runge_rule(method, a, b, eps)
            results.append(result)

        print_results_table(eps, results)


if __name__ == "__main__":
    main()
