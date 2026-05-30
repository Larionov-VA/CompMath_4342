import matplotlib.pyplot as plt

def main():
    with open('res.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    random_n = []
    random_time = []
    hilbert_n = []
    hilbert_rel = []

    for line in lines:
        if line.startswith('2x2') or line.startswith('SINGULAR'):
            continue
        parts = line.split()
        if not parts:
            continue
        if parts[0] == 'R' and len(parts) == 5:
            random_n.append(int(parts[1]))
            random_time.append(float(parts[2]))
        elif parts[0] == 'H':
            hilbert_n.append(int(parts[1]))
            hilbert_rel.append(float(parts[2]))

    if random_n:
        plt.figure()
        plt.plot(random_n, random_time, 'o-')
        plt.xlabel('n')
        plt.ylabel('time (s)')
        plt.title('Время решения случайной матрицы')
        plt.grid(True)
        plt.savefig('time.png')

    if hilbert_n:
        plt.figure()
        plt.semilogy(hilbert_n, hilbert_rel, 'o-')
        plt.xlabel('n')
        plt.ylabel('relative error')
        plt.title('Погрешность матрицы Гильберта')
        plt.grid(True)
        plt.savefig('hilbert.png')
   
    plt.show()

if __name__ == '__main__':
    main()