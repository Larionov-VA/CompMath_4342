from matplotlib import pyplot as plt

with open("./operation_metrix.txt") as file:
    f = file.readlines()

eps, n, root = [], [], []
for line in f:
    deps, dn, droot = line.split()
    eps.append(float(deps))
    n.append(int(dn))
    root.append(float(droot))

plt.figure(figsize=(10, 6))
plt.semilogx(eps, n, marker='o', linestyle='-', color='b')
plt.gca().invert_xaxis()
plt.xlabel('Точность (Eps)')
plt.ylabel('Количество итераций (N)')
plt.title('Зависимость числа итераций от точности')
plt.grid(True, which="both", ls="-")
plt.show()


with open("./root_metrix.txt") as file:
    f2 = file.readlines()

rootM, nM, delta = [], [], []
for line in f2:
    drootM, dnM, ddelta = line.split()
    rootM.append(float(drootM))
    nM.append(int(dnM))
    delta.append(float(ddelta))

plt.figure(figsize=(10, 6))
plt.semilogx(delta, rootM, marker='o', linestyle='-', color='b')
plt.gca().invert_xaxis()
plt.xlabel('Корень (root)')
plt.ylabel('Погрешность округления (delta)')
plt.title('Зависимось значения корня от погрешности округления')
plt.grid(True, which="both", ls="-")
plt.show()