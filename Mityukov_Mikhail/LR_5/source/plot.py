import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


# График 1: число итераций от точности Eps
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
plt.savefig('newton_iterations.png', dpi=220, bbox_inches='tight')
plt.close()

# График 2: погрешность корня от точности округления Delta
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
plt.savefig('newton_sensitivity.png', dpi=220, bbox_inches='tight')
plt.close()