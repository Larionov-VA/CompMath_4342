import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def f(x: float) -> float:
    return math.exp(1.0 / (x * x)) - math.log(x)


def read_csv(path: str):
    with open(path, newline='', encoding='utf-8') as fobj:
        return list(csv.DictReader(fobj))


base = Path('.')
speed_rows = read_csv(base / 'speed_results.csv')
cond_rows = read_csv(base / 'conditioning_results.csv')

# ------------------------------------------------------------
# Рисунок 1: график функции и отрезок отделения корня
# ------------------------------------------------------------
xs = [2.4 + i * (2.0 / 500.0) for i in range(501)]
ys = [f(x) for x in xs]

plt.figure(figsize=(7.2, 4.6))
plt.plot(xs, ys)
plt.axhline(0.0)
plt.axvline(3.0, linestyle='--')
plt.axvline(4.0, linestyle='--')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('f(x) = exp(1/x^2) - ln(x)')
plt.tight_layout()
plt.savefig('figure_function.png', dpi=220)
plt.close()

# ------------------------------------------------------------
# Рисунок 2: число итераций в зависимости от требуемой точности
# ------------------------------------------------------------
xs_speed = [-math.log10(float(row['eps'])) for row in speed_rows]
ys_speed = [int(row['iters']) for row in speed_rows]

plt.figure(figsize=(7.2, 4.6))
plt.plot(xs_speed, ys_speed, marker='o')
plt.xlabel('-log10(Eps)')
plt.ylabel('Число итераций')
plt.title('Метод простых итераций: скорость сходимости')
plt.tight_layout()
plt.savefig('figure_speed.png', dpi=220)
plt.close()

# ------------------------------------------------------------
# Рисунок 3: абсолютная погрешность в зависимости от шага округления
# ------------------------------------------------------------
xs_cond = [float(row['delta']) for row in cond_rows if row['status'] == 'ok']
ys_cond = [float(row['abs_error']) for row in cond_rows if row['status'] == 'ok']

plt.figure(figsize=(7.2, 4.6))
plt.loglog(xs_cond, ys_cond, marker='o')
plt.xlabel('Delta')
plt.ylabel('|x_delta - x_ref|')
plt.title('Чувствительность к округлению в phi(x)')
plt.tight_layout()
plt.savefig('figure_conditioning.png', dpi=220)
plt.close()
