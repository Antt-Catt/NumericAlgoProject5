import numpy as np

def integration_n_trapezoidal(f,a,b,N):
    s=0
    h = (b-a)/N
    for i in range(N):
        s+= (h/2)*(f(a + i*h) + f(a+(i+1)*h))
    return s


def integration_n_midlepoint(f,a,b,N):
    s = 0
    h = (b-a)/N
    for i in range(N):
        s+=h*f(a+i*h + h/2)
    return s

def integration_n_Simson(f,a,b,N):
    return (4/3) *integration_n_trapezoidal(f, a, b, N/2) -  (1/3) * integration_n_trapezoidal(f, a, b, N)