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
