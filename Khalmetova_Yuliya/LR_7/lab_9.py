import math


# 2) Функция для вычисления подынтегральной функции
def f(x):
    return math.sin(x ** 3)

# 1) Программы-функции для вычисления интегралов
def rect_integral(a, b, n):
    """Метод средних прямоугольников"""
    h = (b - a) / n
    total = 0
    for i in range(n):
        x_i = a + i * h
        total += f(x_i + h / 2)
    return h * total


def trap_integral(a, b, n):
    """Метод трапеций"""
    h = (b - a) / n
    total = 0
    for i in range(n):
        x_i = a + i * h
        x_next = a + (i + 1) * h
        total += f(x_i) + f(x_next)
    return (h / 2) * total


def simpson_integral(a, b, n):
    """Метод Симпсона (парабол)"""
    h = (b - a) / n
    m = n // 2  # Количество сдвоенных шагов
    total = 0
    for i in range(m):
        x_2i = a + 2 * i * h
        x_2i_1 = a + (2 * i + 1) * h
        x_2i_2 = a + (2 * i + 2) * h
        total += f(x_2i) + 4 * f(x_2i_1) + f(x_2i_2)
    return (h / 3) * total


# 3) Головная программа с оценкой по Рунге
def solve_with_runge(method_func, method_name, a, b, eps):
    n = 2  # Начинаем с n=2 (четное число, обязательно для Симпсона)
    I_h = method_func(a, b, n)

    while True:
        n *= 2  # Удваиваем число разбиений
        I_h2 = method_func(a, b, n)

        # Разница по модулю
        diff = abs(I_h2 - I_h)

        # Определяем делитель для правила Рунге
        if method_name == "simpson":
            error = diff / 15
        else:
            error = diff / 3

        # Проверяем, достигли ли мы нужной точности
        if error <= eps:
            if method_name == "rect":
                final_val = I_h2 + error
            elif method_name == "trap":
                final_val = I_h2 - error
            else:  # simpson
                final_val = I_h2 - error

            return final_val, n

        # Если точность не достигнута, текущий I_h2 становится старым I_h для следующего шага
        I_h = I_h2


# 4) Проведение вычислений для разных эпсилон
def main():
    a = 1.0
    b = 2.0
    epsilons = [0.01, 0.001, 0.0001,0.00001, 0.000001]

    methods = [
        (rect_integral, "rect", "Прямоугольники"),
        (trap_integral, "trap", "Трапеции"),
        (simpson_integral, "simpson", "Симпсон")
    ]

    for eps in epsilons:
        print(f"\n--- Точность eps = {eps} ---")
        for func, name, label in methods:
            result, final_n = solve_with_runge(func, name, a, b, eps)
            print(f"{label}: Интеграл = {result:.6f}, n = {final_n}")


if __name__ == "__main__":
    main()