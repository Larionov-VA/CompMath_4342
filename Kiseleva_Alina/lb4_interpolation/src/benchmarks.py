import matplotlib.pyplot as plt
import numpy as np
from math import pow

f_spline = open("spline.txt")
f_func = open("func.txt")
f_mid = open("mids.txt")
f_lagr = open("lagrange.txt")

mids_x, mids_y = [], []
arr_x, arr_y = [], []
lagr_vals = [float(s) for s in f_lagr]
spline_coeffs = []

for s in f_spline:
    coeffs = list(map(float, s.split(' ')))
    spline_coeffs.append(coeffs)

for s in f_func:
    x, y = list(map(float, s.split(' ')))
    arr_x.append(x)
    arr_y.append(y)

for s in f_mid:
    x, y = list(map(float, s.split(' ')))
    mids_x.append(x)
    mids_y.append(y)

diff_spline, diff_lagr = [], []
spline_x = []
spline_y = []

for i in range(len(mids_x)):
    spline_val = 0.0
    for j in range(len(spline_coeffs[i])):
        spline_val += spline_coeffs[i][j] * pow(mids_x[i] - arr_x[i], j)

    diff_spline.append(abs(mids_y[i] - spline_val))
    diff_lagr.append(abs(mids_y[i] - lagr_vals[i]))

for i in range(len(arr_x) - 1):
    sub_x = np.linspace(arr_x[i], arr_x[i+1], 10)   #10 точек на отрезке
    for x in sub_x:
        val = 0.0
        for j in range(len(spline_coeffs[i])):
            val += spline_coeffs[i][j] * pow(x - arr_x[i], j)
        spline_x.append(x)
        spline_y.append(val)

f_spline.close()
f_func.close()
f_mid.close()
f_lagr.close()

plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(arr_x, arr_y, 'bo', label='Узлы', markersize=4) # Только точки узлов
plt.plot(spline_x, spline_y, 'r-', label='Кубический сплайн', linewidth=1)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Сравнение функции и сплайна')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
x_fine = np.linspace(min(arr_x), max(arr_x), 200)
y_fine = [np.tan(x) - 1/x for x in x_fine]
plt.plot(x_fine, y_fine, 'g-', label='f(x) = tan(x) - 1/x', linewidth=1.5)
plt.plot(arr_x, arr_y, 'ro', label='Узлы интерполяции', markersize=4)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Исходная функция с узлами')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 3)
plt.semilogy(mids_x, diff_spline, 'b-', marker='o', markersize=3)
plt.xlabel('x (середина отрезка)')
plt.ylabel('|f(x) - S(x)|')
plt.title('Погрешность сплайн-интерполяции')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.semilogy(mids_x, [d + 1e-18 for d in diff_lagr], 'r-', marker='s', markersize=3)
plt.xlabel('x (середина отрезка)')
plt.ylabel('|f(x) - L(x)|')
plt.title('Погрешность интерполяции Лагранжа')
plt.grid(True)

plt.tight_layout()
plt.savefig('plot.png', dpi=150, bbox_inches='tight')