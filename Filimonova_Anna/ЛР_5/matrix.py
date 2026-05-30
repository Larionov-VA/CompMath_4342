import numpy as np

def hilbert(n):
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            H[i, j] = 1.0 / (i + j + 1)
    return H

sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
for n in sizes:
    A = np.random.rand(n, n)
    with open(f'random_{n}.txt', 'w') as f:
        f.write(f'{n}\n')
        for row in A:
            f.write(' '.join(f'{x:.15f}' for x in row) + '\n')
    cond = np.linalg.cond(A)
    print(f'random_{n}.txt cond = {cond:.6e}')

hsizes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
for n in hsizes:
    H = hilbert(n)
    with open(f'hilbert_{n}.txt', 'w') as f:
        f.write(f'{n}\n')
        for row in H:
            f.write(' '.join(f'{x:.15f}' for x in row) + '\n')
    cond = np.linalg.cond(H)
    print(f'hilbert_{n}.txt cond = {cond:.6e}')