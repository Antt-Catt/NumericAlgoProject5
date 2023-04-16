import load_foil as lf
import numpy as np
import matplotlib.pyplot as plt
import part_1 as p1
import part_2 as p2


air_density = 1.293

def f_lambda(l, h, f):
    return lambda x: (1 - l) * f(x) + l * 3 * h

def pression(f,a,b):
    return 1/2 + air_density*p2.L(f,a, b)

if __name__ == "__main__":
    (dim, ex, ey, ix, iy) = lf.load_foil("bacnlf.dat")
    
    ########################################
    # for ex and ey
    n = int(dim[0])
    
    yp1 = 0
    ypn = 0
    
    res = p1.spline(ex, ey, n, yp1, ypn)
    
    xint = np.linspace(ex[0], ex[n - 1], 1000)
    
    fint_supp = p1.interp(ex, ey, res, n)
    
    ########################################
    # for ix and iy
    n = int(dim[1])
    
    res = p1.spline(ix, iy, n, yp1, ypn)
    
    fint_inf = p1.interp(ix, iy, res, n)
    
    ########################################
    hmin = min(iy)
    hmax = max(ey)
    
    x_supp = np.linspace(ex[0], ex[int(dim[0]) - 1], 100)
    y_supp = np.zeros(100)
    
    x_inf = np.linspace(ix[0], ix[int(dim[1]) - 1], 100)
    y_inf = np.zeros(100)
    
    for i in range(20):
        f_supp = f_lambda(i / 20, hmax, fint_supp)
        f_inf = f_lambda(i / 20, hmin, fint_inf)
        
        for j in range(100):
            y_supp[j] = f_supp(x_supp[j])
            y_inf[j] = f_inf(x_inf[j])

        if i == 0:
            plt.plot(x_supp, y_supp)
            plt.plot(x_inf, y_inf)
        else:
            plt.plot(x_supp, y_supp, color = 'k')
            plt.plot(x_inf, y_inf, color = 'k')

    
    
    # plt.ylim(-0.2, 0.2)
    plt.show()
