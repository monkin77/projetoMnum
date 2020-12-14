import math

def rk2(lb, ub, h, pi):
    x = pi[0]
    y = pi[1]
    z = pi[2]
    while(round(x,2) < ub):
        delta1Y = h * dydx(x+h/2, y+(h/2)*dydx(x,y,z), z+(h/2)*dzdx(x,y,z))
        delta1Z = h * dzdx(x+h/2, y+(h/2)*dydx(x,y,z), z+(h/2)*dzdx(x,y,z))
        x += h
        y += delta1Y
        z += delta1Z
        if( (round(x,2) == 0.10) or (round(x, 2) == 0.50) ):
            print(x, y, z)
        #print(i)
    return (x,y,z)

def GaussJacobi(f1, f2, f3, x1, x2, x3, k):     # k = num iteraçoes    
    for _ in range(k):
        xa1 = x1        # Em GaussJacobi utiliza-se o mesmo x,y,z para calcular os proximos valores
        xa2 = x2
        x1 = f1(x2, x3)
        x2 = f2(xa1, x3)
        x3 = f3(xa1, xa2)
    return (x1, x2, x3)

def GaussSeidel(f1, f2, f3, x1, x2, x3, k):
    for _ in range(k):  
        x1 = f1(x2, x3)         # Em GaussSeidel utiliz-se sempre os valores mais atualizados calcular os proximos
        x2 = f2(x1, x3)
        x3 = f3(x1, x2)
    return (x1, x2, x3)

def Trapezio(f, h, x0, xf):
    n = int((xf-x0)//h)
    result = f(x0)+f(xf)
    sum = 0
    for i in range(1, n):
        sum += f(x0 + h*i)
    result = (h/2)*(result + 2*sum)
    return result


def Simpson(f, h, x0, xf):
    n = int((xf-x0)//(2*h))
    sum1 = 0
    sum2 = 0    
    for i in range(1, 2*n, 2):      # not 2n-1 like the formula because on a for loop it goes until 2*n excluding 2*n = i < 2n
        sum1 += f(x0+h*i)
    for i in range(2, 2*n -1, 2):
        sum2 += f(x0 + h*i)
    result = (h/3) * (f(x0) + f(xf) + 4*sum1 +2*sum2)
    return result

def EulerStep(f, x0, xf, y0, h):
    y = y0
    x = x0
    while round(x, 2) < xf:
        y += f(x, y)*h
        x += h
    return y