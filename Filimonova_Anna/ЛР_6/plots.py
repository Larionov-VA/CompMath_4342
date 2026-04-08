import matplotlib.pyplot as plt
import sys

def main():
    if len(sys.argv) < 3:
        print("Использование: plots.py <тип> <выходной_файл>")
        print("тип: eps (итерации от точности) или delta (корень от погрешности)")
        sys.exit(1)

    plot_type = sys.argv[1]
    output = sys.argv[2]

    params = []   
    iters = []    
    roots = []    

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            p = float(parts[0])  
            it = int(parts[1])     
            r = float(parts[2])  
        except ValueError:
            continue
        params.append(p)
        iters.append(it)
        roots.append(r)

    if not params:
        print("Нет данных для построения графика")
        return

    plt.figure(figsize=(8, 6))

    if plot_type == 'eps':
        plt.plot(params, iters, 'o-')
        plt.xscale('log')
        plt.xlabel('eps')
        plt.ylabel('iterations')
    elif plot_type == 'delta':
        plt.plot(params, roots, 'o-')
        plt.xscale('log')
        plt.xlabel('delta')
        plt.ylabel('root')
    else:
        print("Неизвестный тип. Используйте eps или delta")
        return

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(output, dpi=150)
    plt.show()

if __name__ == '__main__':
    main()