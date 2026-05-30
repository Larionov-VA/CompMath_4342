from matplotlib import pyplot as plt

with open("./main_metrix.txt") as file:
    f = file.readlines()

result, eps, n = [], [], []
for line in f:
    droot, deps, dn = line.split()
    eps.append(float(deps))
    n.append(int(dn))
    result.append(float(droot))

plt.figure(figsize=(10, 6))
plt.semilogx(eps, n, marker='o', linestyle='-', color='b')
plt.gca().invert_xaxis()
plt.xlabel('Точность (Eps)')
plt.ylabel('Количество итераций (N)')
plt.title('Зависимость числа итераций от точности')
plt.grid(True, which="both", ls="-")
plt.show()


with open("./conditionality_metrix.txt") as file:
    f2 = file.readlines()

resultM, delta, nM = [], [], []
for line in f2:
    drootM, ddelta, dnM = line.split()
    resultM.append(float(drootM))
    nM.append(int(dnM))
    delta.append(float(ddelta))

plt.figure(figsize=(10, 6))
plt.semilogx(delta, resultM, marker='o', linestyle='-', color='b') 
plt.gca().invert_xaxis() 
plt.xlabel('Погрешность округления (Delta)')
plt.ylabel('Найденный корень')
plt.title('Исследование обусловленности: влияние округления на результат')
plt.grid(True, which="both", ls="-")
plt.show()