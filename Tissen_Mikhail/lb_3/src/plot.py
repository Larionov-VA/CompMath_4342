import csv
import math

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "DejaVu Sans"

LEFT = 3.0
RIGHT = 4.0


def f(x):
    if x <= 0.0:
        raise ValueError("x должен быть больше 0")
    return math.exp(1.0 / (x * x)) - math.log(x)


def bisect(left, right, eps=1e-10):
    f_left = f(left)
    f_right = f(right)

    if f_left * f_right > 0.0:
        raise ValueError("неверный отрезок изоляции")

    while (right - left) >= 2.0 * eps:
        mid = 0.5 * (left + right)
        f_mid = f(mid)

        if f_mid == 0.0:
            return mid
        if f_mid * f_left < 0.0:
            right = mid
        else:
            left = mid
            f_left = f_mid

    return 0.5 * (left + right)


def finish_plot(path):
    plt.savefig(path, dpi=200, bbox_inches="tight")
    if "agg" not in plt.get_backend().lower():
        plt.show()
    plt.close()


root = bisect(LEFT, RIGHT)

# ============================================================
# График 0: сама функция и отрезок изоляции
# ============================================================
xs = [1.0 + i * 0.01 for i in range(401)]
ys = [f(x) for x in xs]

plt.figure(figsize=(9, 5))
plt.plot(xs, ys, color="navy", linewidth=2, label="f(x) = exp(1/x^2) - ln(x)")
plt.axhline(0.0, color="black", linewidth=1)
plt.axvspan(LEFT, RIGHT, color="gold", alpha=0.2, label="Отрезок изоляции [3, 4]")
plt.scatter(
    [LEFT, RIGHT, root],
    [f(LEFT), f(RIGHT), 0.0],
    color=["green", "red", "black"],
    zorder=3,
)
plt.annotate(f"F(3) = {f(LEFT):.4f}", (LEFT, f(LEFT)), xytext=(3.12, 0.22))
plt.annotate(f"F(4) = {f(RIGHT):.4f}", (RIGHT, f(RIGHT)), xytext=(4.08, -0.28))
plt.annotate(f"Корень x* = {root:.6f}", (root, 0.0), xytext=(3.18, -0.52))
plt.xlabel("x", fontsize=11)
plt.ylabel("f(x)", fontsize=11)
plt.title("График функции f(x) = exp(1/x^2) - ln(x)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
finish_plot("function_v18.png")

# ============================================================
# График 1: число итераций от точности Eps
# ============================================================
eps, iters = [], []
with open("iterations_vs_eps.csv") as f_csv:
    next(f_csv)
    for row in csv.reader(f_csv):
        eps.append(float(row[0]))
        iters.append(int(row[1]))

plt.figure(figsize=(8, 5))
plt.semilogx(eps, iters, "bo-", linewidth=2, markersize=6)
plt.xlabel("Точность Eps", fontsize=11)
plt.ylabel("Число итераций N", fontsize=11)
plt.title("Сходимость метода бисекции", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.7)
finish_plot("iterations_v18.png")

# ============================================================
# График 2: погрешность корня от точности округления Delta
# ============================================================
delta, err = [], []
with open("sensitivity_vs_delta.csv") as f_csv:
    next(f_csv)
    for row in csv.reader(f_csv):
        d = float(row[0])
        if d > 0.0:
            delta.append(d)
            err.append(float(row[2]))

plt.figure(figsize=(8, 5))
plt.loglog(delta, err, "ro-", linewidth=2, markersize=6)
plt.xlabel("Шаг округления Delta", fontsize=11)
plt.ylabel("|x - x_опорный|", fontsize=11)
plt.title("Чувствительность к ошибкам округления", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.7)
finish_plot("sensitivity_v18.png")
