import load_foil as lf
import numpy as np
import matplotlib.pyplot as plt
import part_1 as p1
import part_2 as p2

air_density = 1.293

def f_lambda(l, h, f):
    return lambda x: (1 - l) * f(x) + l * 3 * h

def t_lambda(l, h, f):
    return lambda x: -(1 - l) * f(x) - l * 3 * h

def pressure(f, a, b):
    return 1/2 + air_density * (p2.L(f, a, b))**2

if __name__ == "__main__":
    (dim, ex, ey, ix, iy) = lf.load_foil("bacnlf.dat")
    
    ########################################
    # for ex and ey
    n = int(dim[0])
    
    yp1 = 1e30 # to have second derivative = 0
    ypn = 1e30
    
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
    
    x_high = np.linspace(ex[0], ex[int(dim[0]) - 1], 100)
    y_high = np.zeros(100)
    
    x_low = np.linspace(ix[0], ix[int(dim[1]) - 1], 100)
    y_low = np.zeros(100)

    functions_high = []
    functions_low = []
    
    for i in range(20):
        f_high = f_lambda(i / 20, hmax, fint_supp)
        functions_high.append(f_high)
        f_low = f_lambda(i / 20, hmin, fint_inf)
        functions_low.append(f_low)
        
        for j in range(100):
            y_high[j] = f_high(x_high[j])            
            y_low[j] = f_low(x_low[j])
            
        if i == 0:
            plt.plot(x_high, y_high)
            plt.plot(x_low, y_low)
        else:
            plt.plot(x_high, y_high, color = 'k')
            plt.plot(x_low, y_low, color = 'k')

    
    plt.ylim(-0.15, 0.2)
    plt.show()

    press_high = [pressure(f, xint[0], xint[-1]) for f in functions_high]
    press_low = [pressure(f, xint[0], xint[-1]) for f in functions_low]

    X = np.linspace(0, 1, 100)
    Y = np.linspace(-0.15, 0.2, 100)
    Z = np.full((100, 100), min(min(press_high), min(press_low)))

    # Remplissage du tableau 2D avec les coefficients
    for i in range(len(functions_high)):
        for x in X:
            for y in Y:
                if (functions_high[i](x) < y):
                    Z[int((y + 0.15) * (100 - 0) / (0.2 + 0.15)) - 1, int(x * 99)] = press_high[i]

    for i in range(len(functions_low)):
        for x in X:
            for y in Y:
                if (functions_low[i](x) > y):
                    Z[int((y + 0.15) * (100 - 0) / (0.2 + 0.15)), int(x * 99)] = press_low[i]

    # Affichage de la carte des coefficients
    plt.imshow(Z, origin='lower', extent=(0, 1, -0.15, 0.2), aspect='auto', cmap='hot')
    plt.colorbar()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
