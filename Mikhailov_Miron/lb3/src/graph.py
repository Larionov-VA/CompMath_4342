import matplotlib.pyplot as plt

def main():
    # 1. График зависимости числа итераций от точности Eps
    numEps = int(input().strip())
    eps_values = []
    iter_values = []
    roots_eps = []

    for _ in range(numEps):
        parts = input().split()
        eps = float(parts[0])
        n = float(parts[1])
        x = float(parts[2])
        eps_values.append(eps)
        iter_values.append(n)
        roots_eps.append(x)

    # Построение графика 1
    plt.figure()
    plt.plot(eps_values, iter_values, marker='o', linestyle='-', color='b')
    plt.xscale('log')
    plt.xlabel('Точность Eps')
    plt.ylabel('Число итераций')
    plt.title('Зависимость числа итераций от точности Eps')
    plt.grid(True)
    plt.savefig('eps_n.png', dpi=150)
    plt.close()

    # Вывод в консоль
    print("\nИсследование зависимости числа операций от точности eps")
    for i in range(numEps):
        print(f"eps = {eps_values[i]:g}, n = {iter_values[i]:g}, x = {roots_eps[i]:.12f}")

    # 2. График зависимости корня от погрешности Delta (чувствительность)
    numDelta = int(input().strip())
    delta_values = []
    iter_delta = []
    root_values = []

    for _ in range(numDelta):
        parts = input().split()
        delta = float(parts[0])
        n = float(parts[1])
        x = float(parts[2])
        delta_values.append(delta)
        iter_delta.append(n)
        root_values.append(x)

    # Построение графика 2
    plt.figure()
    plt.plot(delta_values, root_values, marker='s', linestyle='-', color='r')
    plt.xscale('log')               # Delta меняется на порядки – логарифмическая шкала
    plt.xlabel('Погрешность Delta')
    plt.ylabel('Приближённый корень')
    plt.title('Зависимость вычисленного корня от погрешности Delta')
    plt.grid(True)
    plt.savefig('delta_x.png', dpi=150)
    plt.close()

    # Вывод в консоль
    print("\nИсследование чувствительности к ошибкам исходных данных")
    for i in range(numDelta):
        print(f"delta = {delta_values[i]:g}, n = {iter_delta[i]:g}, x = {root_values[i]:.12f}")

if __name__ == "__main__":
    main()