from matplotlib import pyplot as plt

def plot(x_data, y_data, x_label, y_label, name, file_name, is_log_xscale = False):
    plt.plot(x_data, y_data)
    plt.grid()
    if is_log_xscale: 
        plt.xscale('log')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(name)
    plt.savefig(file_name)
    plt.show()

def main():
    eps_values = []
    n_values = []
    delta_values = []
    x_values = []

    print("Исследование зависимости числа операций от точности eps");

    numEps = int(input())
    for _ in range(numEps):
        eps, n, x = [float(el) for el in input().split()]
        print(f"eps = {eps}, n = {n}, x = {x}")
        eps_values.append(eps)
        n_values.append(n)

    plot(eps_values, n_values, "Точность Eps", "Количество итераций", "Зависимость числа итераций от точности Eps", "eps_n.png", is_log_xscale=True)

    print("\nИсследование чувствительности к ошибкам исходных данных")
    numDelta = int(input())
    for _ in range(numDelta):
        delta, n, x = [float(el) for el in input().split()]
        print(f"delta = {delta}, n = {n}, x = {x}")
        delta_values.append(delta)
        x_values.append(x)

    plot(delta_values, x_values, "Погрешность Delta", "Приближенный корень", "Зависимость вычисленного корня от погрешности Delta", "delta_x.png")

if __name__ == "__main__":
    main()