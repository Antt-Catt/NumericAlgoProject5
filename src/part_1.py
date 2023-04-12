import load_foil as lf
import numpy as np
import matplotlib.pyplot as plt

# This function is the translation of the function in the Numerical Recipies in python
# Input :
# x and y such as f(x[i]) = y[i] with f the function to interpolate
# yp1 the first derivative of the unterpolating function at point 1
# ypn the first derivative of the unterpolating function at point n
def spline(x, y, n, yp1, ypn):
    u = np.zeros(n)

    y2 = np.zeros(n)

    if yp1 > 0.99e30:
        y2[0] = u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)

    for i in range(1, n-1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p

    if ypn > 0.99e30:
        qn = un = 0.0
    else:
        qn = 0.5
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]))

    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1)

    for k in range(n-2, -1, -1):
        y2[k] = y2[k] * y2[k+1] + u[k]

    return y2


def splint(xa, ya, y2a, n, x):
    def nrerror(error_text):
        raise ValueError(error_text)

    res = np.zeros(n)

    klo = 0
    khi = n-1
    
    while (khi - klo > 1):
        k = (khi + klo) // 2
        if (xa[k] > x):
            khi = k
        else:
            klo = k
    
    h = xa[khi] - xa[klo]
    if (h == 0):
        nrerror("Bad xa input to routine splint")
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h

    y = a * ya[klo] + b * ya[khi] + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * (h ** 2) / 6
    return y

if __name__ == "__main__":

    def f(x):
        return x**3

    def g(x):
        return x**5

    def h(x):
        return f(x) + g(x)

    def k(x):
        return x**4 + 3*x**5

    def interp(xa, ya, y2, n):
        return lambda x: splint(xa, ya, y2, n, x)

    n = 50
    x = np.linspace(-1, 1, n)


    ##################################################
    # tests for f(x)
    y = np.array(list(map(f, x)))

    yp1 = 3
    ypn = 3

    res = spline(x, y, n, yp1, ypn)

    fint = interp(x, y, res, n)
    xint = np.linspace(-1, 1, 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
        if abs(yint[i] - f(xint[i])) > 1e-5:
            print("Not enough accurate")

    plt.plot(xint, yint, label = "f(x) = x**3")

    ##################################################
    # tests for g(x)
    y = np.array(list(map(g, x)))

    yp1 = 5
    ypn = 5

    res = spline(x, y, n, yp1, ypn)

    fint = interp(x, y, res, n)
    xint = np.linspace(-1, 1, 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
        if abs(yint[i] - g(xint[i])) > 1e-5:
            print("Not enough accurate")

    
    plt.plot(xint, yint, label = "g(x) = x**5")


    ##################################################
    # tests for h(x)
    y = np.array(list(map(h, x)))

    yp1 = 8
    ypn = 8

    res = spline(x, y, n, yp1, ypn)

    fint = interp(x, y, res, n)
    xint = np.linspace(-1, 1, 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
        if abs(yint[i] - h(xint[i])) > 1e-5:
            print("Not enough accurate")

    plt.plot(xint, yint, label = "h(x) = f(x) + g(x)")


    ##################################################
    # tests for k(x)
    y = np.array(list(map(k, x)))

    yp1 = 11
    ypn = 19

    res = spline(x, y, n, yp1, ypn)

    fint = interp(x, y, res, n)
    xint = np.linspace(-1, 1, 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
        if abs(yint[i] - k(xint[i])) > 1e-5:
            print("Not enough accurate")

    plt.plot(xint, yint, label = "k(x) = x**4 + 3*x**5")

    plt.legend()
    plt.show()


    ##################################################
    # tests for airfoil
    (dim, ex, ey, ix, iy) = lf.load_foil("bacnlf.dat")

    plt.scatter(ex, ey)
    plt.scatter(ix, iy)
    plt.ylim(-0.3, 0.3)

    ########################################
    # for ex and ey
    n = int(dim[0])
    
    x = ex
    y = ey
    
    yp1 = (y[1] - y[0]) / (x[1] - x[0])
    ypn = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])

    res = spline(x, y, n, yp1, ypn)
        
    fint = interp(x, y, res, n)
    xint = np.linspace(x[0], x[n - 1], 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
    
    plt.plot(xint, yint)
    print(yint)

    ########################################
    # for ix and iy
    n = int(dim[1])
    
    x = ix
    y = iy

    yp1 = (y[1] - y[0]) / (x[1] - x[0])
    ypn = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    
    res = spline(x, y, n, yp1, ypn)
    
    fint = interp(x, y, res, n)
    xint = np.linspace(x[0], x[n - 1], 1000)
    yint = np.zeros(1000)
    
    for i in range(1000):
        yint[i] = fint(xint[i])
    
    plt.plot(xint, yint)
    
    plt.show()
