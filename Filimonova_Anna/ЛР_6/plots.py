import matplotlib.pyplot as plt
import math

eps_neg_log = []
iterations = []

with open('result.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) == 2:
            try:
                e = float(parts[0])
                it = int(parts[1])
                if it >= 0:
                    eps_neg_log.append(-math.log10(e))
                    iterations.append(it)
            except ValueError:
                continue

plt.figure(figsize=(8, 5))
plt.plot(eps_neg_log, iterations, 'o-', linewidth=2, markersize=8)
plt.xlabel('$-\log_{10} \epsilon$')
plt.ylabel('Число итераций')
plt.title('Зависимость числа итераций от точности')
plt.grid(True)
plt.savefig('iter_eps.png', dpi=150)
plt.show()