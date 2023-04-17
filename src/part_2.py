import numpy as np
import matplotlib.pyplot as plt
import math

def integration_n_trapezoidal(f,a,b,N):
    s=(f(a) +f(b))*0.5
    h = (b-a)/N
    for i in range(1,N):
        s+= f(a + i*h)
    return h*s


def integration_n_midlepoint(f,a,b,N):
    s = 0
    h = (b-a)/N
    for i in range(N):
        s+=h*f(a+i*h + h/2)
    return s

def integration_n_Simpson(f,a,b,N):
    return (4/3) *integration_n_trapezoidal(f, a, b, N//2) -  (1/3) * integration_n_trapezoidal(f, a, b, N)

def integration_n_left(f, a, b, N):
    h = (b-a)/N
    integral = 0
    for i in range(N):
        integral += f(a + i*h)
    integral *= h
    return integral

def integration_n_right(f, a, b, N):
    h = (b-a)/N
    integral = 0
    for i in range(1, N+1):
        integral += f(a + i*h)
    integral *= h
    return integral

def integration_epsilon(f,a,b,eps,meth):
    n = 2
    Nmax = 1000
    while n < Nmax:
        if (abs(meth(f,a,b,2*n)- meth(f,a,b,n)) < eps):
            return meth(f,a,b,n)
        n *= 2
    return None


def polynomial(x):
    return x**3 - 3*x + 2

a = 0
b = 2

exact_integral = 2

# Test trapezoidal rule
N = 10
approx_integral = integration_n_trapezoidal(polynomial, a, b, N)
error = abs(approx_integral - exact_integral)
print(f"Trapezoidal rule with N={N}: Approximate integral = {approx_integral}, Error = {error}")

# Test midpoint rule
N = 10
approx_integral = integration_n_midlepoint(polynomial, a, b, N)
error = abs(approx_integral - exact_integral)
print(f"Midpoint rule with N={N}: Approximate integral = {approx_integral}, Error = {error}")

# Test Simpson's rule
N = 10
approx_integral = integration_n_Simpson(polynomial, a, b, N)
error = abs(approx_integral - exact_integral)
print(f"Simpson's rule with N={N}: Approximate integral = {approx_integral}, Error = {error}")

# Test left endpoint rule
N = 10
approx_integral = integration_n_left(polynomial, a, b, N)
error = abs(approx_integral - exact_integral)
print(f"Left endpoint rule with N={N}: Approximate integral = {approx_integral}, Error = {error}")

# Test right endpoint rule
N = 10
approx_integral = integration_n_right(polynomial, a, b, N)
error = abs(approx_integral - exact_integral)
print(f"Right endpoint rule with N={N}: Approximate integral = {approx_integral}, Error = {error}")

# Test adaptive integration
eps = 0.01
approx_integral = integration_epsilon(polynomial, a, b, eps, integration_n_Simpson)
error = abs(approx_integral - exact_integral)
print(f"Adaptive integration with epsilon={eps}: Approximate integral = {approx_integral}, Error = {error}")



# Convergence graphs

X = [ k  for k in range(2,50)]
Exact = [2 for i in range(2,50)]
Ytrap = [ integration_n_trapezoidal(polynomial, a, b, x) for x in X]
Ymid = [integration_n_midlepoint(polynomial, a, b, x) for x in X]
Ysimp = [integration_n_Simpson(polynomial, a, b, x) for x in X]
Yleft = [integration_n_left(polynomial, a, b, x) for x in X]
Yright = [integration_n_right(polynomial, a, b, x) for x in X]

plt.plot(X,Ytrap,label = "Trapezoidal")
plt.plot(X,Ymid,label = "Midpoint")
plt.plot(X, Ysimp, label = "Simpson")
plt.plot(X,Yleft,label = "Right rectangle")
plt.plot(X,Yright, label = "Left rectangle")
plt.plot(X,Exact,label="exact integral", color ="black")


plt.legend()
plt.ylim([1.5,2.5])
plt.show()

# Length of plane curves 
def L(f,a,b,method=integration_n_Simpson):
    def derivative(g,x):
        h = 1e-3
        return (g(x+h) - g(x))/(h) 
    def function(x):
        return math.sqrt(1+derivative(f, x)**2)
    return method(function,a,b,10)

    

