import pandas as pd
import matplotlib.pyplot as plt

# 1) Iterations vs Eps
eps_df = pd.read_csv("iterations_vs_eps_newton.csv", sep=";")

plt.figure()
plt.plot(eps_df["Eps"], eps_df["Iterations"], marker="o", linewidth=1)
plt.xscale("log")
plt.xlabel("Eps")
plt.ylabel("Iterations")
plt.title("Newton: Iterations vs Eps")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("iterations_vs_eps_newton.png", dpi=200)
plt.show()

# 2) Error vs Delta (Delta=0 убираем для log по X)
delta_df = pd.read_csv("sensitivity_vs_delta_newton.csv", sep=";")
delta_df_nonzero = delta_df[delta_df["Delta"] > 0].copy()

plt.figure()
plt.plot(delta_df_nonzero["Delta"], delta_df_nonzero["AbsErrorToRef"], marker="o", linewidth=1)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Delta")
plt.ylabel("|x_delta - x_ref|")
plt.title("Newton: Sensitivity (Error vs Delta)")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("sensitivity_vs_delta_newton.png", dpi=200)
plt.show()