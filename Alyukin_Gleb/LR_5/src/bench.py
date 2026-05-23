import matplotlib.pyplot as plt
import csv
import numpy as np

# 1. Функция f(x) и ее производная f'(x)
def f(x):
    return x**2 - 5.0 * np.sin(x)

def df(x):
    return 2.0 * x - 5.0 * np.cos(x)

# 2. Метод Ньютона для подсчета итераций
def newton_count(x0, eps):
    n = 0
    x = x0
    while True:
        y = f(x)
        y1 = df(x)
        if y1 == 0: # Защита от деления на ноль
            return n
        dx = y / y1
        x = x - dx
        n += 1
        # Условие выхода из цикла как в твоем C++ коде: while (fabs(DX) > Eps)
        if abs(dx) <= eps:
            break
    return n

# 3. Подготовка данных в порядке от 10^-7 до 10^-1
# (от самой высокой точности к самой низкой)
eps_values = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
data = []

x0 = 3.0 # Выбрано по условию f(x0)*f''(x0) > 0

for e in eps_values:
    n_iters = newton_count(x0, e)
    data.append([f"{e:.7f}", n_iters])

# 4. Сохранение в CSV файл
csv_filename = 'lab5_results.csv'
with open(csv_filename, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Eps', 'Iterations'])
    writer.writerows(data)

print(f"Данные сохранены в {csv_filename}")

# 5. Построение графика
x_eps = [float(row[0]) for row in data]
y_iters = [row[1] for row in data]

plt.figure(figsize=(10, 6))
# Строим график метода Ньютона (синим цветом)
plt.plot(x_eps, y_iters, marker='o', linestyle='-', color='blue', label='Метод Ньютона')

# Настройка осей
plt.xscale('log') # Логарифмическая шкала для Eps
plt.grid(True, which="both", linestyle=':', alpha=0.7)

# Оформление
plt.title('Зависимость количества итераций от величины Eps\n(Метод Ньютона, x0 = 3.0)', fontsize=14)
plt.xlabel('Eps (логарифмическая шкала: от 10^-7 до 10^-1)', fontsize=12)
plt.ylabel('Количество итераций (N)', fontsize=12)

# Добавляем текстовые подписи над точками
for i in range(len(x_eps)):
    plt.annotate(str(y_iters[i]), (x_eps[i], y_iters[i]),
                 textcoords="offset points", xytext=(0,10), ha='center')

plt.legend()
plt.show()