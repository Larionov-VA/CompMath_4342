import random
import math
import time
import numpy as np
import matplotlib.pyplot as plt


# Блок 1: генераторы матриц
def generate_random_matrix(n, interval=100):
    """Случайная матрица"""
    return [[random.uniform(-interval, interval) for _ in range(n)] for _ in range(n)]  # создание вложенных списков n x n со случайными числами

def generate_hilbert_matrix(n):
    """Матрица Гильберта (плохая)"""
    return [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)]  # генерация элементов по формуле Гильберта для анализа погрешностей

def generate_diagonal_dominant_matrix(n, interval=100):
    """Диагонально-доминантная матрица (идеальная)"""
    A = generate_random_matrix(n, interval)
    for i in range(n):
        row_sum = sum(abs(A[i][j]) for j in range(n) if i != j)  # сумма модулей всех элементов строки кроме диагонального
        A[i][i] = row_sum + random.uniform(1.0, 10.0)           # делаем диагональный элемент больше суммы остальных
        if random.choice([True, False]):
            A[i][i] = -A[i][i]                                  # рандомная смена знака диагонали для общности
    return A

# Блок 2: подготовка системы и метод Гаусса
def prepare_system(A):
    """Генерирует секретный ответ X_exact и вычисляет столбец B"""
    n = len(A)
    X_exact = [random.uniform(-10, 10) for _ in range(n)]       # создаем эталонный вектор ответов для проверки точности
    B = [0.0] * n
    for i in range(n):
        B[i] = sum(A[i][j] * X_exact[j] for j in range(n))      # вычисляем правую часть уравнения b через известный X
    return B, X_exact

def gaussian_elimination(A_input, B_input):
    """Метод Гаусса с частичным выбором по столбцу"""
    n = len(A_input)
    A = [row[:] for row in A_input]  # копирование матрицы, чтобы не изменять оригинал
    B = B_input[:]                   # копирование вектора b

    # Прямой ход
    for k in range(n):
        max_row = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[max_row][k]):
                max_row = i  # поиск максимального элемента в столбце для обеспечения численной устойчивости

        if max_row != k:
            A[k], A[max_row] = A[max_row], A[k]  # перестановка строк в матрице A
            B[k], B[max_row] = B[max_row], B[k]  # соответствующая перестановка в векторе b

        if abs(A[k][k]) < 1e-12:
            continue  # пропуск итерации, если ведущий элемент близок к нулю (защита от вырожденности)

        pivot = A[k][k]
        for j in range(k, n):
            A[k][j] /= pivot  # нормировка ведущей строки (деление на диагональный элемент)
        B[k] /= pivot

        for i in range(k + 1, n):
            factor = A[i][k]
            for j in range(k, n):
                A[i][j] -= factor * A[k][j]  # обнуление элементов под диагональю
            B[i] -= factor * B[k]            # применение тех же операций к вектору b

    # Обратный ход
    X_calc = [0.0] * n
    for i in range(n - 1, -1, -1):
        X_calc[i] = B[i]
        for j in range(i + 1, n):
            X_calc[i] -= A[i][j] * X_calc[j]  # последовательное нахождение неизвестных от последнего к первому

    return X_calc

def calculate_error(X_exact, X_calc):
    """Вычисляет относительную погрешность (невязку)"""
    n = len(X_exact)
    diff_norm = math.sqrt(sum((X_exact[i] - X_calc[i]) ** 2 for i in range(n)))  # Евклидова норма разности векторов
    exact_norm = math.sqrt(sum(x ** 2 for x in X_exact))                       # Евклидова норма эталонного вектора
    return diff_norm / (exact_norm + 1e-15)                                    # расчет относительной погрешности (невязки)


# Блок 3: режим исследования (таблица и график)
def run_research():
    sizes = [10, 20, 50, 100, 200]
    matrix_types = [
        ("Random", "Случайная", generate_random_matrix, 'blue'),
        ("Diag", "Диагонально-доминантная", generate_diagonal_dominant_matrix, 'green'),
        ("Hilbert", "Гильберта", generate_hilbert_matrix, 'red')
    ]

    results = {
        "Random": {"N": [], "Time": [], "Error": [], "Cond": []},
        "Diag": {"N": [], "Time": [], "Error": [], "Cond": []},
        "Hilbert": {"N": [], "Time": [], "Error": [], "Cond": []}
    }

    print("\nЗапуск вычислений... Это может занять пару секунд.\n")

    print(f"{'N':<6} | {'Тип матрицы':<25} | {'Время (мс)':<12} | {'Относ. невязка':<20} | {'Число обусловленности'}")
    print("-" * 95)

    for n in sizes:
        for key, name, generator, color in matrix_types:
            if key == "Hilbert" and n > 50:
                continue  # пропускаем матрицы Гильберта больших порядков из-за колоссальной погрешности

            A = generator(n)
            B, X_exact = prepare_system(A)

            start_time = time.perf_counter()
            X_calc = gaussian_elimination(A, B)
            elapsed_time = (time.perf_counter() - start_time) * 1000  # замер времени решения в миллисекундах

            error = calculate_error(X_exact, X_calc)
            cond = np.linalg.cond(np.array(A))  # расчет числа обусловленности матрицы средствами numpy

            results[key]["N"].append(n)
            results[key]["Time"].append(elapsed_time)
            results[key]["Error"].append(error)
            results[key]["Cond"].append(cond)

            print(f"{n:<6} | {name:<25} | {elapsed_time:<12.4f} | {error:<20.5e} | {cond:.5e}")
        print("-" * 95)

    # Отрисовка графиков
    plt.style.use('seaborn-v0_8-whitegrid')

    # 1. Время
    plt.figure(figsize=(10, 6))
    plt.plot(results["Random"]["N"], results["Random"]["Time"], marker='o', color='purple', linewidth=2)
    plt.title('Зависимость времени выполнения от размера матрицы N', fontsize=14)
    plt.xlabel('Размер матрицы (N)', fontsize=12)
    plt.ylabel('Время выполнения (мс)', fontsize=12)
    plt.grid(True)
    plt.savefig('Time_Graph.png', dpi=300, bbox_inches='tight')  # сохранение графика времени выполнения

    # 2. Невязка
    plt.figure(figsize=(10, 6))
    for key, name, _, color in matrix_types:
        if results[key]["N"]:
            plt.plot(results[key]["N"], results[key]["Error"], marker='s', color=color, label=name, linewidth=2)
    plt.yscale('log')  # использование логарифмической шкалы для наглядности погрешности
    plt.title('Зависимость относительной невязки от размера матрицы N', fontsize=14)
    plt.xlabel('Размер матрицы (N)', fontsize=12)
    plt.ylabel('Относительная погрешность', fontsize=12)
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig('Error_Graph.png', dpi=300, bbox_inches='tight')  # сохранение графика невязки

    # 3. Обусловленность
    plt.figure(figsize=(10, 6))
    for key, name, _, color in matrix_types:
        if results[key]["N"]:
            plt.plot(results[key]["N"], results[key]["Cond"], marker='^', color=color, label=name, linewidth=2)
    plt.yscale('log')
    plt.title('Зависимость числа обусловленности от размера матрицы N', fontsize=14)
    plt.xlabel('Размер матрицы (N)', fontsize=12)
    plt.ylabel('Число обусловленности cond(A)', fontsize=12)
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig('Condition_Graph.png', dpi=300, bbox_inches='tight')  # сохранение графика числа обусловленности

    print("\nУСПЕХ! Графики сохранены в твоей папке проекта: Time_Graph.png, Error_Graph.png, Condition_Graph.png")

# Блок 4: решение конкретной системы
def run_solve():
    print("\n--- РЕШЕНИЕ КОНКРЕТНОЙ СИСТЕМЫ ---")
    try:
        n = int(input("Введите размер матрицы N: "))  # пользовательский ввод размерности
        if n < 2:
            print("Матрица должна быть хотя бы 2x2!")
            return

        print("\nВыберите тип матрицы:")
        print("1 - Случайная")
        print("2 - Диагонально-доминантная")
        print("3 - Матрица Гильберта")
        m_type = int(input("Ваш выбор (1/2/3): "))

        if m_type == 1:
            name = "Случайная"
            A = generate_random_matrix(n)
        elif m_type == 2:
            name = "Диагонально-доминантная"
            A = generate_diagonal_dominant_matrix(n)
        elif m_type == 3:
            name = "Матрица Гильберта"
            A = generate_hilbert_matrix(n)
        else:
            print("Ошибка: неверный тип матрицы.")
            return

        B, X_exact = prepare_system(A)

        start_time = time.perf_counter()
        X_calc = gaussian_elimination(A, B)
        elapsed_time = (time.perf_counter() - start_time) * 1000

        error = calculate_error(X_exact, X_calc)
        cond = np.linalg.cond(np.array(A))

        print(f"\n[РЕЗУЛЬТАТЫ] Матрица: {name} | Размер: {n}x{n}")
        print(f"Время Гаусса: {elapsed_time:.4f} мс")
        print(f"Относительная погрешность: {error:.5e}")
        print(f"Число обусловленности:     {cond:.5e}\n")

        print(f"{'Индекс':<8} | {'Эталонный X (Секретный)':<25} | {'Мой ответ (Гаусс)'}")
        print("-" * 65)

        limit = min(n, 10)  # ограничение вывода, чтобы не забивать экран при больших N
        for i in range(limit):
            print(f"x[{i + 1:<2}]  | {X_exact[i]:<25.5f} | {X_calc[i]:.5f}")

        if n > 10:
            print(f"... и еще {n - 10} элементов скрыто.")

    except ValueError:
        print("Ошибка: Пожалуйста, вводите только целые числа.")

# Блок 5: главное меню (терминал)
def main():
    while True:
        print("\n" + "=" * 45)
        print(" ЛАБОРАТОРНАЯ РАБОТА: МЕТОД ГАУССА")
        print("=" * 45)
        print("1. Запустить исследование")
        print("2. Решить конкретную систему (Вывод x1, x2...)")
        print("0. Выйти из программы")

        try:
            choice = input("\nВаш выбор: ")
            if choice == '1':
                run_research()  # вызов автоматического исследования
            elif choice == '2':
                run_solve()     # вызов ручного решения СЛАУ
            elif choice == '0':
                print("Программа завершена.")
                break
            else:
                print("Неизвестная команда. Введите 1, 2 или 0.")
        except KeyboardInterrupt:
            print("\nПрограмма принудительно завершена.")
            break


if __name__ == "__main__":
    main()  # точка входа в программу