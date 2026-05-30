import random
import time

# Блок 1: Генераторы матриц

def make_dominant_matrix(dim):
    """Генерирует правильную (диагонально-доминантную) матрицу"""
    mat = [[0.0] * dim for _ in range(dim)]
    for i in range(dim):
        row_sum = 0.0
        for j in range(dim):
            if i != j:
                val = random.uniform(-10.0, 10.0)
                mat[i][j] = val
                row_sum += abs(val)

        mat[i][i] = row_sum + random.uniform(1.0, 5.0)
        if random.choice([True, False]):
            mat[i][i] = -mat[i][i]

    return mat


def make_random_matrix(dim):
    """Генерирует абсолютно случайную матрицу"""
    mat = [[0.0] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            mat[i][j] = random.uniform(-10.0, 10.0)
        if mat[i][i] == 0:
            mat[i][i] = 1.0
    return mat


def make_vectors(mat, dim):
    """Создает секретный ответ и правую часть уравнения"""
    exact_x = [random.uniform(-10.0, 10.0) for _ in range(dim)]
    rhs_b = [0.0] * dim

    for i in range(dim):
        temp_sum = 0.0
        for j in range(dim):
            temp_sum += mat[i][j] * exact_x[j]
        rhs_b[i] = temp_sum

    return exact_x, rhs_b

# Блок 2: Метод Якоби

def solve_by_jacobi(dim, mat, rhs, exact_x):
    """Возвращает: статус, число итераций, финальную ошибку"""
    x_current = [0.0] * dim
    x_next = [0.0] * dim
    iters = 0

    while True:
        x_last = x_current[:]

        try:
            for i in range(dim):
                s = rhs[i]
                for j in range(dim):
                    if i != j:
                        s -= mat[i][j] * x_last[j]
                x_next[i] = s / mat[i][i]

            max_err = 0.0
            for k in range(dim):
                err = abs(x_next[k] - x_last[k])
                if err > max_err:
                    max_err = err

            x_current[:] = x_next[:]
            iters += 1

            if max_err > 1e20:
                return "РАСХОДИТСЯ", iters, max_err

            if iters >= 10000:
                return "ПРЕВЫШЕН ЛИМИТ", iters, max_err

            if max_err <= 1e-3:
                return "УСПЕХ", iters, max_err

        except OverflowError:
            return "ФАТАЛЬНОЕ РАСХОЖДЕНИЕ", iters, float('inf')

# Блок 3: Тестирование и Сравнение

def print_matrix(mat, rhs, dim, name):
    print(f"\n{name}:")
    print("-" * 65)
    for i in range(dim):
        row_str = "  ".join([f"{mat[i][j]:8.2f}" for j in range(dim)])
        print(f"| {row_str} |  =  {rhs[i]:8.2f}")
    print("-" * 65)

# Блок 4: Новое исследование

def calculate_norm_B(mat, dim):
    """Считает норму ||B|| для Таблицы 1"""
    max_norm = 0.0
    for i in range(dim):
        row_sum = 0.0
        for j in range(dim):
            if i != j:
                row_sum += abs(mat[i][j])
        # Норма по строке = сумма модулей недиагональных элементов / диагональный элемент
        current_norm = row_sum / abs(mat[i][i])
        if current_norm > max_norm:
            max_norm = current_norm
    return max_norm


def run_research():
    """Проводит тесты и записывает их в файл research_results.txt"""
    sizes = [10, 50, 100, 200, 500, 1000, 10000]
    filename = "research_results.txt"

    with open(filename, "w", encoding="utf-8") as f:
        f.write("Исследование обусловленности\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'N':<5} | {'Итерации':<10} | {'Время (сек)':<12} | {'Погрешность':<12} | {'||B||':<8}\n")
        f.write("-" * 80 + "\n")

        print("Запуск исследования... Это может занять несколько секунд.")
        for n in sizes:
            mat = make_dominant_matrix(n)
            exact, b = make_vectors(mat, n)
            norm_b = calculate_norm_B(mat, n)

            start = time.perf_counter()
            status, iters, final_err = solve_by_jacobi(n, mat, b, exact)
            end = time.perf_counter()

            line = f"{n:<5} | {iters:<10} | {end - start:<12.5f} | {final_err:<12.2e} | {norm_b:.4f}\n"
            f.write(line)
            print(f"Размер {n} готов.")

    print(f"\nГотово! Результаты сохранены в файл: {filename}")


#main с выбором режима
def main():
    print("Выберите режим:")
    print("1. Решить одну систему уравнений")
    print("2. Исследование для отчета (Запись в файл)")

    choice = input("Ваш выбор: ")

    if choice == "2":
        run_research()
        return
    try:
        size = int(input("Введите размер системы: "))
    except ValueError:
        print("Ошибка ввода.")
        return

    print("\n" + "=" * 50)
    print(" ЭКСПЕРИМЕНТ 1: ДИАГОНАЛЬНО-ДОМИНАНТНАЯ МАТРИЦА")
    print("=" * 50)
    mat_good = make_dominant_matrix(size)
    exact_good, b_good = make_vectors(mat_good, size)
    print_matrix(mat_good, b_good, size, "Диагонально-доминантная матрица:")

    status, iters, err = solve_by_jacobi(size, mat_good, b_good, exact_good)
    print(f"Статус:     {status}")
    print(f"Итераций:   {iters}")
    print(f"Погрешность: {err:.5e}")

    print("\n" + "=" * 50)
    print(" ЭКСПЕРИМЕНТ 2: СЛУЧАЙНАЯ МАТРИЦА")
    print("=" * 50)
    mat_bad = make_random_matrix(size)
    exact_bad, b_bad = make_vectors(mat_bad, size)
    print_matrix(mat_bad, b_bad, size, "Матрица составленная из случайных значений:")

    status, iters, err = solve_by_jacobi(size, mat_bad, b_bad, exact_bad)
    print(f"Статус:     {status}")
    print(f"Итераций:   {iters}")
    print(f"Погрешность: {err:.5e}")


if __name__ == "__main__":
    main()