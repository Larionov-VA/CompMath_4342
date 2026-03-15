import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Функция варианта 18
# ============================================================
def f(x: np.ndarray) -> np.ndarray:
    return np.exp(1.0 / (x * x)) - np.log(x)

# ============================================================
# График 0: функция и интервал изоляции корня
# ============================================================
x = np.linspace(2.5, 4.2, 600)
y = f(x)

plt.figure(figsize=(8, 5))
plt.plot(x, y)
plt.axhline(0.0)
plt.axvline(3.0, linestyle="--")
plt.axvline(4.0, linestyle="--")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Вариант 18: f(x) = exp(1/x^2) - ln(x)")
plt.grid(True)
plt.tight_layout()
plt.savefig("chords_isolation_plot.png", dpi=200)
plt.close()

# ============================================================
# График 1: число итераций от точности Eps
# ============================================================
eps_df = pd.read_csv("eps_iter_profile.csv", sep=";").sort_values("Eps")
plt.figure(figsize=(8, 5))
plt.plot(eps_df["Eps"], eps_df["Iterations"], marker="o", linewidth=1)
plt.xscale("log")
plt.xlabel("Eps")
plt.ylabel("Iterations")
plt.title("Метод хорд: итерации vs Eps")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("chords_eps_curve.png", dpi=200)
plt.close()

# ============================================================
# График 2: погрешность корня от точности округления Delta
# ============================================================
delta_df = pd.read_csv("delta_impact_profile.csv", sep=";")
delta_ok = delta_df[(delta_df["Ok"] == 1) & (delta_df["Delta"] > 0)].copy()
delta_ok = delta_ok.sort_values("Delta")

plt.figure(figsize=(8, 5))
plt.plot(delta_ok["Delta"], delta_ok["AbsErrorToRef"], marker="o", linewidth=1)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Delta")
plt.ylabel("|x_delta - x_ref|")
plt.title("Метод хорд: чувствительность к округлению")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("chords_delta_curve.png", dpi=200)
plt.close()
