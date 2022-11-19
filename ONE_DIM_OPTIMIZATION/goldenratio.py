import numpy as np
from math import sqrt
import matplotlib.pyplot as plt


def f(x: float) -> float:
    return x ** 2 + 3 * x * (np.log(x) - 1)


def golden_ratio(a: float, b: float, eps: float) -> tuple:
    k: int = 1
    phi: float = (1 + sqrt(5)) / 2
    x1: float = b - (b - a) / phi
    x2: float = a + (b - a) / phi
    while abs(b - a) > eps:
        y1: float = f(x1)
        y2: float = f(x2)
        if y1 >= y2:
            a = x1
            x1 = x2
            x2 = a + (b - a) / phi
        else:
            b = x2
            x2 = x1
            x1 = b - (b - a) / phi
        k += 1
    xmin: float = (a + b) / 2
    return xmin, k


if __name__ == '__main__':
    a: float = 0.5
    b: float = 1.0
    eps: float = 5e-5
    print(f'xmin = {golden_ratio(a, b, eps)[0]}')
    print(f'k = {golden_ratio(a, b, eps)[1]}')
    x: np.ndarray = np.linspace(1e-3, 3, 200)
    y: np.ndarray = f(x)
    
    plt.title('Golden Ratio')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.plot(x, y, label=f'f(x)\nx_min = {golden_ratio(a, b, eps)[0]}\ny_min = {f(golden_ratio(a, b, eps)[0])}\nk = {golden_ratio(a, b, eps)[1]}')
    
    xmin: float = golden_ratio(a, b, eps)[0]
    ymin: float = f(xmin)

    plt.scatter(xmin, ymin)
    plt.text(xmin, ymin - 0.5, 'min', weight='bold')

    plt.grid()
    plt.legend()
    plt.show()