import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 24,          
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'legend.fontsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'figure.figsize': (24, 18) 
})

data = np.loadtxt("plot_data.dat")
x = data[:, 0]
f = data[:, 1]
poly = data[:, 2]
spl = data[:, 3]

plt.figure()
plt.plot(x, f, 'k-', linewidth=3, label='f(x)')
plt.plot(x, poly, 'r--', linewidth=2, label='Полином Ньютона')
plt.plot(x, spl, 'b-.', linewidth=2, label='Кубический сплайн')
f_min, f_max = np.min(f), np.max(f)
margin = (f_max - f_min) * 0.3
plt.ylim(f_min - margin, f_max + margin)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Сравнение интерполяции')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plot.png")
plt.show()