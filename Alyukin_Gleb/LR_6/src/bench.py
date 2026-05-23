import matplotlib.pyplot as plt
import csv
import numpy as np

# 1. Преобразованная функция phi(x) = sqrt(5 * sin(x))
def phi(x):
    # Проверка области определения
    val = 5.0 * np.sin(x)
    if val < 0:
        return 0
    return np.sqrt(val)

# 2. Метод простых итераций для подсчета количества шагов
def iter_count(x0, eps):
    n = 0
    x_old = x0
    max_limit = 1000 # Защита от бесконечного цикла
    
    while n < max_limit:
        x_new = phi(x_old)
        n += 1
        # Условие выхода: разность между соседними приближениями меньше Eps
        if abs(x_new - x_old) < eps:
            break
        x_old = x_new
        
    return n

# 3. Подготовка данных (от высокой точности 10^-7 к низкой 10^-1)
eps_values = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
data = []

# Начальное приближение x0 = 2.1 (выбрано на этапе отделения корней)
x0 = 2.1 

for e in eps_values:
    n_iters = iter_count(x0, e)
    data.append([f"{e:.7f}", n_iters])

# 4. Сохранение результатов в CSV файл
csv_filename = 'lab6_results.csv'
with open(csv_filename, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Eps', 'Iterations'])
    writer.writerows(data)

print(f"Данные сохранены в {csv_filename}")

# 5. Построение графика
x_eps = [float(row[0]) for row in data]
y_iters = [row[1] for row in data]

plt.figure(figsize=(10, 6))
# Строим график метода простых итераций (зеленым цветом)
plt.plot(x_eps, y_iters, marker='d', linestyle='-', color='green', label='Метод простых итераций')

# Настройка осей
plt.xscale('log') # Логарифмическая шкала для Eps
plt.grid(True, which="both", linestyle=':', alpha=0.7)

# Оформление
plt.title('Зависимость количества итераций от величины Eps\n(Метод простых итераций, x0 = 2.1, q ≈ 0.62)', fontsize=14)
plt.xlabel('Eps (логарифмическая шкала: от 10^-7 до 10^-1)', fontsize=12)
plt.ylabel('Количество итераций (N)', fontsize=12)

# Добавляем текстовые подписи над точками
for i in range(len(x_eps)):
    plt.annotate(str(y_iters[i]), (x_eps[i], y_iters[i]),
                 textcoords="offset points", xytext=(0,10), ha='center')

plt.legend()
plt.show()