
import matplotlib.pyplot as plt
import sys

N, t_e, e_e, t_g, e_g = [], [], [], [], []
try:
    with open('log.txt', 'r') as f:
        for line in f:
            data = list(map(float, line.split()))
            N.append(data[0]); t_e.append(data[1]); e_e.append(data[2]); t_g.append(data[3]); e_g.append(data[4])
except Exception as e:
    print('Ошибка чтения:', e)
    sys.exit(1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
ax1.plot(N, t_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax1.plot(N, t_g, 'r-s', markersize=4, label='GaussianElimination')
ax1.set_yscale('log'); ax1.grid(True, linestyle='--'); ax1.set_title('Время выполнения'); ax1.set_xlabel('N'); ax1.set_ylabel('Время, мс'); ax1.legend()

ax2.plot(N, e_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax2.plot(N, e_g, 'r-s', markersize=4, label='GaussianElimination')
ax2.set_yscale('log'); ax2.grid(True, linestyle='--'); ax2.set_title('Относительная погрешность'); ax2.set_xlabel('N'); ax2.set_ylabel('Погрешность'); ax2.legend()

plt.tight_layout()
plt.savefig('graphs.png', dpi=300)
plt.show()
