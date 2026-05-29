import matplotlib
# matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import sys

data = open(sys.argv[1], "r").readlines()
xs = []
comps = []
precs = []

for line in data:
    if line:
        x, y1, y2 = list(map(float, line.split()))
        xs.append(x)
        comps.append(y1)
        precs.append(y2 + 2)


plt.plot(xs, comps, "bo", label="Time (ms)")
plt.plot(xs, comps, "k--")
plt.title("Complexity ~ Table size (N)")
plt.show()

plt.plot(xs, precs, "bo", label="Precision in decimal places")
plt.plot(xs, precs, "k--")
plt.title("Precision ~ Table size (N)")
plt.show()


