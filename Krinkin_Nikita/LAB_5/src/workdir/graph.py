import matplotlib
matplotlib.use("TkAgg")

import os
import sys
from matplotlib import pyplot as plt


if __name__ == "__main__":
    data_x = []
    data_y = []


    if not os.path.exists(sys.argv[1]):
        raise FileNotFoundError(f"No {sys.argv[1]} file!")
    file_nm = sys.argv[1]
    file = open(file_nm, "r")

    if sys.argv[2] == "precision":
        for line in file.readlines():
            if line:
                x, y = line.split()
                data_x.append(int(x))
                data_y.append(int(y))

        plt.plot(data_x, data_y, "bo", label="Result")
        plt.plot(data_x, data_y, "k--")

        plt.plot(data_x, list(map(lambda x: 2 * x**(1/2), data_x)), "r--", label="Abstract f(x) = sqrt(x)")
        plt.title("Iterations ~ Precision")
        plt.xlabel("Precision, number of decimal places")
        plt.ylabel("Number of iterations")

    else:
        for line in file.readlines():
            if line:
                x, y = line.split()
                data_x.append(float(x))
                data_y.append(int(y))

        plt.axvline(x=0.70711, color='red', linestyle='--', label='f\'(x) = 0 \nx = max(f\'(x))')
        plt.plot(data_x, data_y, "bo", label="Result")
        # plt.plot(data_x, data_y, "k--")

        plt.title("Iterations ~ first approximation")
        plt.xlabel("Value of x0 -- first approximation")
        plt.ylabel("Number of iterations")
    
    plt.legend()
    plt.show()