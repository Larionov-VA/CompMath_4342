import pandas as pd
import matplotlib.pyplot as plt

# 1) Итерации от Eps
eps_df = pd.read_csv("iterations_vs_eps.csv")

plt.figure()
plt.plot(eps_df["Eps"], eps_df["Iterations"], marker="o", linewidth=1)
plt.xscale("log")
plt.xlabel("Eps")
plt.ylabel("Iterations")
plt.title("Bisection: Iterations vs Eps")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("iterations_vs_eps.png", dpi=200)
plt.show()

# 2) Чувствительность: ошибка от Delta (для Delta=0 лог-шкала не подходит)
delta_df = pd.read_csv("sensitivity_vs_delta.csv")

# Уберём Delta = 0, чтобы можно было поставить log по X
delta_df_nonzero = delta_df[delta_df["Delta"] > 0].copy()

plt.figure()
plt.plot(delta_df_nonzero["Delta"], delta_df_nonzero["AbsErrorToRef"], marker="o", linewidth=1)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Delta")
plt.ylabel("|x - x_ref|")
plt.title("Sensitivity: Error vs Delta (Eps fixed)")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("sensitivity_vs_delta.png", dpi=200)
plt.show()