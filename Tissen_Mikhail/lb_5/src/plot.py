import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

LEFT = 3.0
RIGHT = 4.0
X0 = 3.0


# ============================================================
# Функция варианта 18
# ============================================================
def f(x: float) -> float:
    if x <= 0.0:
        raise ValueError("x must be positive")
    return math.exp(1.0 / (x * x)) - math.log(x)


# ============================================================
# Первая производная функции варианта 18
# ============================================================
def df(x: float) -> float:
    if x <= 0.0:
        raise ValueError("x must be positive")
    return -2.0 * math.exp(1.0 / (x * x)) / (x ** 3) - 1.0 / x


# ============================================================
# Один шаг метода Ньютона
# ============================================================
def newton_step(x: float) -> float:
    return x - f(x) / df(x)


# ============================================================
# График 0: функция и геометрический смысл метода Ньютона
# ============================================================
xs = [2.2 + i * 0.005 for i in range(451)]
ys = [f(x) for x in xs]
x1 = newton_step(X0)
x2 = newton_step(x1)
root = newton_step(x2)

plt.figure(figsize=(8.2, 5.4))
plt.plot(xs, ys, color='navy', linewidth=2.2, label='f(x)', zorder=2)
plt.axhline(0.0, color='black', linewidth=1.0, zorder=1)
plt.axvspan(LEFT, RIGHT, alpha=0.09, label='отрезок [3;4]')

# Точки графика
plt.plot([LEFT, RIGHT], [f(LEFT), f(RIGHT)], 'o', markersize=8)
plt.plot([X0], [f(X0)], 'o', markersize=8)
plt.plot([x1, x2, root], [0.0, 0.0, 0.0], 'o', markersize=6)

# Касательная в x0
k0 = df(X0)
y0 = f(X0)
tx0 = [2.88, x1 + 0.025]
ty0 = [y0 + k0 * (tx0[0] - X0), y0 + k0 * (tx0[1] - X0)]
plt.plot(
    tx0,
    ty0,
    '--',
    color='crimson',
    linewidth=2.2,
    dashes=(8, 4),
    alpha=0.9,
    zorder=5,
    path_effects=[pe.Stroke(linewidth=3.6, foreground='white'), pe.Normal()],
    label='касательная в x₀',
)

# Касательная в x1
k1 = df(x1)
y1 = f(x1)
tx1 = [x1 - 0.06, x2 + 0.05]
ty1 = [y1 + k1 * (tx1[0] - x1), y1 + k1 * (tx1[1] - x1)]
plt.plot(
    tx1,
    ty1,
    ':',
    color='darkorange',
    linewidth=2.4,
    alpha=0.9,
    zorder=5,
    path_effects=[pe.Stroke(linewidth=3.8, foreground='white'), pe.Normal()],
    label='касательная в x₁',
)

plt.scatter(
    [X0, x1],
    [y0, y1],
    s=58,
    color=['crimson', 'darkorange'],
    edgecolors='white',
    linewidths=0.9,
    zorder=6,
)

# Вертикали
plt.vlines(X0, 0.0, y0, colors='gray', linestyles='dotted', linewidth=0.9, alpha=0.9, zorder=1)
plt.vlines(x1, 0.0, y1, colors='gray', linestyles='dotted', linewidth=0.9, alpha=0.9, zorder=1)

# Подписи
label_box = dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.9)
plt.annotate('f(3)', (LEFT, f(LEFT)), xytext=(-22, 12), textcoords='offset points', bbox=label_box, zorder=7)
plt.annotate('f(4)', (RIGHT, f(RIGHT)), xytext=(-12, -18), textcoords='offset points', bbox=label_box, zorder=7)
plt.annotate('x₀ = 3', (X0, 0.0), xytext=(-4, -30), textcoords='offset points', ha='center', bbox=label_box, zorder=7)
plt.annotate('x₁', (x1, 0.0), xytext=(-10, -14), textcoords='offset points', ha='center', bbox=label_box, zorder=7)
plt.annotate('x₂', (x2, 0.0), xytext=(10, -31), textcoords='offset points', ha='center', bbox=label_box, zorder=7)
plt.annotate('x*', (root, 0.0), xytext=(26, -8), textcoords='offset points', ha='center', bbox=label_box, zorder=7)

plt.xlim(2.65, 4.15)
plt.ylim(-0.40, 0.12)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Функция f(x) и первые шаги метода Ньютона')
plt.grid(True, linestyle='--', alpha=0.55)
plt.legend(loc='lower left')
plt.savefig('newton_function_v18.png', dpi=220, bbox_inches='tight')
plt.close()

# ============================================================
# График 1: число итераций от точности Eps
# ============================================================
rows = []
with open('newton_iterations_vs_eps.csv', encoding='utf-8') as fcsv:
    reader = csv.DictReader(fcsv)
    for row in reader:
        rows.append((float(row['Eps']), int(row['Iterations'])))
rows.sort(key=lambda t: t[0])
eps = [r[0] for r in rows]
iters = [r[1] for r in rows]

plt.figure(figsize=(8.0, 5.0))
plt.semilogx(eps, iters, marker='o', linewidth=2.0)
plt.xlabel('Точность Eps')
plt.ylabel('Число итераций N')
plt.title('Сходимость метода Ньютона')
plt.grid(True, which='both', linestyle='--', alpha=0.55)
plt.savefig('newton_iterations_v18.png', dpi=220, bbox_inches='tight')
plt.close()

# ============================================================
# График 2: погрешность корня от точности округления Delta
# ============================================================
rows = []
with open('newton_sensitivity_vs_delta.csv', encoding='utf-8') as fcsv:
    reader = csv.DictReader(fcsv)
    for row in reader:
        d = float(row['Delta'])
        if d > 0.0:
            rows.append((d, float(row['Error'])))
rows.sort(key=lambda t: t[0])
delta = [r[0] for r in rows]
err = [r[1] for r in rows]

plt.figure(figsize=(8.0, 5.0))
plt.loglog(delta, err, marker='o', linewidth=2.0)
plt.xlabel('Точность округления Delta')
plt.ylabel('|x - x_ref|')
plt.title('Чувствительность метода Ньютона к округлению')
plt.grid(True, which='both', linestyle='--', alpha=0.55)
plt.savefig('newton_sensitivity_v18.png', dpi=220, bbox_inches='tight')
plt.close()
