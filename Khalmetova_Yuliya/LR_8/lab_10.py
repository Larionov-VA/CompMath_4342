import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange, CubicSpline
import sys

def run_lab(n_nodes, func_type='runge'):
    """
    Универсальная функция для интерполяции.
    Параметр func_type принимает значения:
      - 'runge': Функция Рунге на отрезке [-1, 1]
      - 'custom': Функция arctg(x) - 1 / x на отрезке [1, 2]
    """

    # 1. НАСТРОЙКА ПАРАМЕТРОВ ПО ФЛАГУ
    if func_type == 'runge':
        a, b = -1.0, 1.0

        def target_func(x):
            return 1 / (1 + 25 * x ** 2)

        func_label = 'Рунге: 1 / (1 + 25x^2)'
        y_limits = (-1.5, 2.5)

    elif func_type == 'custom':
        a, b = 1.0, 2.0

        def target_func(x):
            return np.arctan(x) - 1 / x

        func_label = 'arctg(x) - 1/x'
        y_limits = None

    else:
        print("Ошибка: неизвестный тип функции!")
        return

    print(f"\n{'=' * 50}")
    print(f"РАСЧЕТ: Функция {func_label} | Узлов: N = {n_nodes}")
    print(f"{'=' * 50}")

    # 2. ПОДГОТОВКА СЕТКИ
    x_nodes = np.linspace(a, b, n_nodes)
    y_nodes = target_func(x_nodes)

    x_dense = np.linspace(a, b, 1000)
    y_dense = target_func(x_dense)

    # 3. ИНТЕРПОЛЯЦИЯ ЛАГРАНЖА
    poly_lagrange = lagrange(x_nodes, y_nodes)
    y_lagrange = poly_lagrange(x_dense)

    print("\n[Аналитический аппарат Лагранжа]")
    print(f"Степень: {poly_lagrange.order}")
    print("Коэффициенты (от старшей степени к младшей):")
    print(np.round(poly_lagrange.coef, 4))

    # 4. КУБИЧЕСКИЙ СПЛАЙН
    cs = CubicSpline(x_nodes, y_nodes, bc_type='natural')
    y_spline = cs(x_dense)

    print("\n[Аналитический аппарат Кубического Сплайна]")
    print(f"Количество отрезков: {len(x_nodes) - 1}")
    for i in range(min(3, len(x_nodes) - 1)):
        coefs = np.round(cs.c[:, i], 4)
        print(f"Отрезок [{x_nodes[i]:.2f}, {x_nodes[i + 1]:.2f}]: "
              f"{coefs[0]}x^3 + {coefs[1]}x^2 + {coefs[2]}x + {coefs[3]}")
    if len(x_nodes) - 1 > 3:
        print("... (остальные скрыты)")

    # БЛОК: СРАВНИТЕЛЬНОЕ ИНТЕГРИРОВАНИЕ
    print("\n[СРАВНЕНИЕ РЕЗУЛЬТАТОВ ИНТЕГРИРОВАНИЯ]")

    # 1. Интеграл Лагранжа
    poly_integ = poly_lagrange.integ()
    lagrange_integral = poly_integ(b) - poly_integ(a)
    print(f"Интеграл (через полином Лагранжа):  {lagrange_integral:.6f}")

    # 2. Интеграл Сплайна (ручной метод)
    manual_spline_integral = 0.0
    for i in range(len(x_nodes) - 1):
        h = x_nodes[i + 1] - x_nodes[i]
        a_coef, b_coef, c_coef, d_coef = cs.c[:, i]
        area_i = (a_coef / 4) * h ** 4 + (b_coef / 3) * h ** 3 + (c_coef / 2) * h ** 2 + d_coef * h
        manual_spline_integral += area_i
    print(f"Интеграл (через ручной сплайн):     {manual_spline_integral:.6f}")

    # 3. Интеграл Сплайна (scipy)
    scipy_spline_integral = cs.integrate(a, b)
    print(f"Интеграл (через CubicSpline.integ): {scipy_spline_integral:.6f}")

    # Вывод разницы (для наглядности)
    diff = abs(lagrange_integral - scipy_spline_integral)
    print(f"Разница методов:                   {diff:.6f}")
    # --------------------------------------------------------

    # 5. ОТРИСОВКА ГРАФИКА
    plt.figure(figsize=(9, 5))
    plt.plot(x_dense, y_dense, 'k-', linewidth=4, alpha=0.2, label=f'Оригинал: {func_label}')
    plt.plot(x_dense, y_lagrange, 'r--', linewidth=2, label='Лагранж')
    plt.plot(x_dense, y_spline, 'b-.', linewidth=2, label='Кубический сплайн')
    plt.plot(x_nodes, y_nodes, 'ko', markersize=6, label='Узлы')

    plt.title(f'Сравнение интерполяций (Узлов: {n_nodes}) | Режим: {func_type}', fontsize=12)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend()

    if y_limits:
        plt.ylim(y_limits)

    plt.show()

# ПЕРЕНАПРАВЛЕНИЕ ВЫВОДА В ФАЙЛ
with open('output_results.txt', 'w', encoding='utf-8') as f:
    sys.stdout = f  # Теперь всё, что идет в print(), запишется в файл

    # ЗАПУСК ТЕСТОВ
    print("\n>>> ТЕСТИРУЕМ ФУНКЦИЮ РУНГЕ <<<")
    run_lab(n_nodes=7, func_type='runge')
    run_lab(n_nodes=15, func_type='runge')
    run_lab(n_nodes=21, func_type='runge')
    run_lab(n_nodes=28, func_type='runge')

    print("\n>>> ТЕСТИРУЕМ СВОЮ ФУНКЦИЮ <<<")
    # run_lab(n_nodes=2, func_type='custom')
    # run_lab(n_nodes=4, func_type='custom')
    # run_lab(n_nodes=10, func_type='custom')
    # run_lab(n_nodes=12, func_type='custom')
    # run_lab(n_nodes=100, func_type='custom')

# Возвращаем стандартный вывод обратно в консоль после завершения
sys.stdout = sys.__stdout__
print("Все текстовые результаты успешно сохранены в файл 'output_results.txt'")