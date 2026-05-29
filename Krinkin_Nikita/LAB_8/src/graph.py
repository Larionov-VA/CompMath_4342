import matplotlib
matplotlib.use("TkAgg")   
from matplotlib import pyplot as plt
import sys

pres = open(sys.argv[1], "r").readlines()
comp = open(sys.argv[2], "r").readlines()

xs_p = []
xs_c = []
ys_p = []
ys_c = []

for line1, line2 in zip(pres, comp):
    if line1:
        xp, yp = list(map(int, line1.split()))
        xs_p.append(xp)
        ys_p.append(yp)
    if line2:
        xc, yc = list(map(int, line2.split()))
        xs_c.append(xc)
        ys_c.append(yc)

plt.plot(xs_p, ys_p, "bo", label="Result")
plt.plot(xs_p, ys_p, "k--")
plt.title("Iterations ~ Precision")
plt.xlabel("Precision, number of decimal places")
plt.ylabel("Number of iterations")
plt.legend()
plt.show()

plt.plot(xs_c, ys_c, "bo", label="Result")
plt.plot(xs_c, ys_c, "k--")
plt.title("Iterations ~ Size")
plt.xlabel("Precision, number of decimal places")
plt.ylabel("Number of iterations")
plt.legend()
plt.show()
