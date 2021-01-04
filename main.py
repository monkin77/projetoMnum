import math
import matplotlib.pyplot as plt
import numpy as np

tmax = 1.5 * 60  # unidade : min
doseDiaria = 250    # unidade: mg
doseComp = 250 / 2  # unidade: mg
freq = 12 * 60          # comprimido de 125mg de 12 em 12h (12*60 min)
duracaoTratamento = 7 * 24 * 60  # unidade: minutos
Vap = 3650  # unidade: ml   (Volume de plasma)
Ket = 0.2 / 60   # unidade: m^(-1)   (Taxa de eliminaçao)


def calcAbsorcao(Ka):
    return Ka * ((math.e)**(-Ka * tmax)) - Ket * ((math.e)**(-Ket * tmax))


def calcAbsorcaol(Ka):
    return -tmax * Ka * ((math.e)**(-Ka * tmax))

def g(x):           # Para método Picard-peano
    return Ket * ( (math.e)**(tmax * (x - Ket)) )

# ALGORITMOS DE CALCULAR ZEROS DE FUNÇOES

def newton(f, fl, stopCondition, guess):
    nextX = guess
    while(True):
        #print(nextX, f(nextX))
        x = nextX   # update the value of x (previous x)
        nextX = x - f(x) / fl(x)    # update the value of nextX

        if(abs((nextX - x) / nextX) <= stopCondition):
            return nextX

def bissecao(f, lb, ub, epsilon):
    while( abs( (lb-ub)/lb ) >= epsilon ):
        xMed = (lb+ub)/2
        med = f( xMed)
        if( f(lb) * med < 0 ):
            ub = xMed
        elif( f(lb) * med > 0):
            lb = xMed
        else:
            return xMed
    return lb

def picardPeano(g, epsilon, guess):
    x = guess
    while(True):
        if( abs(g(x) - x) < epsilon ):
            return x
        x = g(x)

Ka1 = bissecao(calcAbsorcao, 0.015, 0.035, 0.0001)
Ka = newton(calcAbsorcao, calcAbsorcaol, 0.0001, 0.01)         # Ka da-me um valor inválido a partir do método de Newton
Ka2 = picardPeano(g, 0.00001, 0.026)
print("Ka1:", Ka1, "Ka2:", Ka, "Ka3:", Ka2)

# 2nd STEP
Mi = 0      # massa no compartimento central
Mp = 0      # massa no compartimento plasmático


def D(t):
    iterationT = t % freq
    if(iterationT <= tmax):     # dose está a ser administrada
        return 3.086 * (10**(-2)) * iterationT
    else:
        return 0

def D2(t):
    iterationT = t %  freq
    if(iterationT <= tmax):     # dose está a ser administrada
        return 1.54 * (10**(-2)) * iterationT
    elif(iterationT <= 2*tmax):
        return -1.54 * (10**(-2)) * iterationT + 2.77
    else:
        return 0

halfLife = 210  #min
Dmax = doseComp / 345    # calculated in the notebook
linearSlope = Dmax / tmax   #slope of the linear part of D(t)
print("Dmax: ", Dmax, "linearSlope: ", linearSlope)

def D_Tooth(t):
    if t <= tmax:
        return linearSlope * t
    elif t >= tmax:
        return Dmax * ( math.e**(-Ket * (t - tmax) ) )   # (t - tmax since the decrase only starts after reaching the max)

def D3(t):
    if t < 0: 
        return 0
    elif 0 <= t <= freq:
        return D_Tooth(t)
    else:
        totalSum = 0
        counter = 1
        while(t > 0):
            if counter > (duracaoTratamento / freq):   # 14 - Number of times the pill is taken
                break
            counter += 1
            totalSum += D_Tooth(t)
            t -= freq
        return totalSum


def dmidt(t, mi):
    return D3(t) - Ka * mi


def dmpdt(t, mi, mp):
    return Ka * mi - Ket * mp

# ALGORITMOS DE RESOLVER SISTEMAS DE EQ. DIFERENCIAIS

# This functions finds the right h to use (Qc <= 2)
def findEulerH(lb, ub, h, pi, dydx, dzdx):    
    while(True):
        s = euler(lb, duracaoTratamento, h, pi, dmidt, dmpdt)
        sl = euler(lb, duracaoTratamento, h/2, pi, dmidt, dmpdt)
        sll = euler(lb, duracaoTratamento, h/4, pi, dmidt, dmpdt)
        y = s[1]
        yl = sl[1]
        yll = sll[1]
        z = s[2]
        zl = sl[2]
        zll = sll[2]
        if( ( (yl - y) / (yll - yl) <= 2 ) and ( (zl - z) / (zll - z) <= 2 ) ):    # Found a good h value -> can calculate error -> if quoc. convergencia <= 2
            absErrY = yll - yl
            absErrZ = zll - zl
            return (absErrY, absErrZ, h)
        h /= 2

def euler(lb, ub, h, pi, dydx, dzdx):
    x = pi[0]   # tempo
    y = pi[1]   # massa no compartimento central
    z = pi[2]   # massa no comp. plasmático
    it = 0
    while(round(x, 2) < ub):
        miArray.append(y)
        mpArray.append(z)
        deltaY = h * dydx(x, y)
        deltaZ = h * dzdx(x, y, z)
        x += h
        y += deltaY
        z += deltaZ
        it += 1
    print("Euler iterations:", it)
    return (x, y, z)

def improvedEuler(lb, ub, h, pi, dydx, dzdx):
    prevPoint = (pi[0], pi[1], pi[2])
    miArray.append(prevPoint[1])
    mpArray.append(prevPoint[2])
    currPoint = [prevPoint[0] + h, prevPoint[1] + dydx(prevPoint[0], prevPoint[1]) * h, prevPoint[2] + dzdx(prevPoint[0], prevPoint[1], prevPoint[2])*h]
    it = 0
    while(round(currPoint[0], 2) < ub):
        miArray.append(currPoint[1])
        mpArray.append(currPoint[2])
        currdydx = dydx(currPoint[0], currPoint[1])
        currdzdx = dydx(currPoint[0], currPoint[1])
        nextPoint = (prevPoint[0] + 2*h, prevPoint[1] + currdydx * 2*h, prevPoint[2] + currdzdx * 2*h )
        futuredydx = dydx(nextPoint[0], nextPoint[1])
        futuredzdx = dzdx(nextPoint[0], nextPoint[1], nextPoint[2])
        deltaY = h * (futuredydx + currdydx) / 2
        deltaZ = h * (futuredzdx + currdzdx) / 2
        prevPoint = currPoint
        currPoint[0] += h
        currPoint[1] += deltaY
        currPoint[2] += deltaZ
        it += 1
    print("Improved Euler iterations:", it)
    return (currPoint[0], currPoint[1], currPoint[2])

def rk2(lb, ub, h, pi, dydx, dzdx):
    x, y, z = pi
    it = 0
    while round(x, 2) < ub:
        miArray.append(y)
        mpArray.append(z)
        delta1Y = h * dydx(x + h / 2, y + (h / 2) * dydx(x, y))
        delta1Z = h * dzdx(x + h / 2, y + (h / 2) * dydx(x, y), z + (h / 2) * dzdx(x, y, z))
        x += h
        y += delta1Y
        z += delta1Z
        it += 1
    print("RK2 iterations:", it)
    return x, y, z

def rk4(lb, ub, h, pi, dydx, dzdx):
    x = pi[0]
    y = pi[1]
    z = pi[2]
    it = 0
    while(round(x, 2) < ub):
        miArray.append(y)
        mpArray.append(z)
        delta1Y = h * dydx(x, y)
        delta1Z = h * dzdx(x, y, z)
        delta2Y = h * dydx(x+h/2, y+delta1Y/2)
        delta2Z = h * dzdx(x+h/2, y+delta1Y/2, z+delta1Z/2)
        delta3Y = h * dydx(x + h/2, y + delta2Y / 2)
        delta3Z = h * dzdx(x + h/2, y + delta2Y / 2, z + delta2Z / 2)
        delta4Y = h * dydx(x+h, y + delta3Y)
        delta4Z = h * dzdx(x+h, y + delta3Y, z + delta3Z)
        deltaY = delta1Y/6 + delta2Y/3 + delta3Y/3 + delta4Y/6
        deltaZ = delta1Z/6 + delta2Z/3 + delta3Z/3 + delta4Z/6
        x += h
        y += deltaY
        z += deltaZ
        it += 1
    print("RK4 iterations:", it)
    return (x, y, z)


lb = 0
ub = duracaoTratamento
# x = tempo / y = Massa no compartimento central / z = Massa no compartimento plasmático
pi = (0, 0, 0)
h = 1
miArray = []
mpArray = []

x = np.arange(0, duracaoTratamento, h)     # range(0, duracaoTratamento, h)

#print(findEulerH(lb, duracaoTratamento, h, pi, dmidt, dmpdt))  -> h to use with Euler method is 1

"""          #D(t) em funçao do tempo
x = np.arange(0, duracaoTratamento, 1)     # range(0, 90, 1)
y = [ D3(i) for i in x ]   # list comprehension
plt.plot(x,y)
plt.xlabel('time / min')
plt.ylabel('D(t) / mg min-1')

# giving a title to my graph
plt.title('Função de Administração')
plt.show()    """

# MASSA DE AMOXICILINA NO PLASMA  &  MASSA DE AMOXICILINA NO COMPARTIMENTO CENTRAL
# y = [ rk4(lb, i , h, pi, dmidt, dmpdt)[2] for i in x ]   # list comprehension
rk4(lb, duracaoTratamento, h, pi, dmidt, dmpdt)       # Much faster to call the function once and append all y to an array
#improvedEuler(lb, duracaoTratamento, h, pi, dmidt, dmpdt)

# CONCENTRAÇAO DE AMOXICILINA NO PLASMA
cmpl = list(map(lambda x: x/Vap, mpArray))
#print("Concentraçao no plasma", cmi)

# CONCENTRAÇAO DE AMOXICILINA NO COMP. CENTRAL
#cmi = list(map(lambda x: x/Vap, miArray))
#print("Concentraçao no plasma", cmi)

# 3rd STEP - Plotting the functions

# plotting the points
plt.plot(x, cmpl)

# naming the x axis
plt.xlabel('concentraçao (mg / ml)')
# naming the y axis
plt.ylabel('massa (mg) ')

# giving a title to my graph
plt.title('Concentraçao plasmatica de amoxicilina')

# function to show the plot
plt.show() 

