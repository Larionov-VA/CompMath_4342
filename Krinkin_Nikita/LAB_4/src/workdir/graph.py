import matplotlib
matplotlib.use("TkAgg")

import os
import sys
from matplotlib import pyplot as plt

data_x = []
data_y = []

if not os.path.exists(sys.argv[1]):
    raise FileNotFoundError(f"No {sys.argv[1]} file!")
file_nm = sys.argv[1]
file = open(file_nm, "r")

for line in file.readlines():
    if line:
        x, y = line.split()
        data_x.append(int(x))
        data_y.append(int(y))

plt.plot(data_x, data_y, "bo", label="Result")
plt.plot(data_x, data_y, "k--")
plt.plot(data_x, list(map(lambda x: 2*x, data_x)), "r--", label="Abstract f(x) = kx")
plt.title("Iterations ~ Precision")
plt.xlabel("Precision, number of decimal places")
plt.ylabel("Number of iterations")
plt.legend()
plt.show()
