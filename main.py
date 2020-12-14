import math
import matplotlib.pyplot as plt 
import numpy as np

tmax = 1.5 * 60  # unidade : min
doseDiaria = 250    # unidade: mg
doseComp = 250 / 2  #unidade: mg
freq = 12 * 60          # comprimido de 125mg de 12 em 12h (12*60 min)
duracaoTratamento = 7 * 24 * 60  # unidade: minutos
Vap = 3650  # unidade: ml   (Volume de plasma)
Ket = 0.2 / 60   # unidade: m^(-1)   (Taxa de eliminaçao) 

def calcAbsorcao(Ka):
    return Ka * ( (math.e)**(-Ka * tmax) ) - Ket* ( (math.e)**(-Ket * tmax) )

def calcAbsorcaol(Ka):
    return -tmax * Ka * ( (math.e)**(-Ka * tmax) )

def newton(f, fl, stopCondition, guess):        
    nextX = guess
    while(True):
        #print(nextX, f(nextX))
        x = nextX   # update the value of x (previous x)
        nextX = x - f(x) / fl(x)    # update the value of nextX

        if( abs( (nextX - x) / nextX ) <= stopCondition):
            return nextX


Ka = newton(calcAbsorcao, calcAbsorcaol, 0.0001, 0.01)
print("Ka:", Ka)

# 2nd STEP
Mi = 0      # massa no compartimento central
Mp = 0      # massa no compartimento plasmático

def D(t):
    iterationT = t % freq
    if( iterationT <= tmax ):     # dose está a ser administrada
        return 3.086 * (10**(-2)) * iterationT
    else:
        return 0
    
def dmidt(t, mi):
    return D(t) - Ka * mi

def dmpdt(t, mi, mp):
    return Ka * mi - Ket * mp

def rk4(lb,ub,h,pi, dydx, dzdx):
    x = pi[0]
    y = pi[1]
    z = pi[2]
    while(round(x, 2) < ub):
        miArray.append(y)
        mpArray.append(z)
        delta1Y = h * dydx(x, y)
        delta1Z = h * dzdx(x, y, z)
        delta2Y = h * dydx(x+h/2, y+delta1Y/2)
        delta2Z = h * dzdx(x+h/2, y+delta1Y/2, z+delta1Z/2 )
        delta3Y = h * dydx(x + h/2, y + delta2Y / 2)
        delta3Z = h * dzdx(x+ h/2, y + delta2Y / 2, z + delta2Z / 2)
        delta4Y = h * dydx(x+h, y + delta3Y)
        delta4Z = h * dzdx(x+h, y + delta3Y, z + delta3Z)
        deltaY = delta1Y/6 + delta2Y/3 + delta3Y/3 + delta4Y/6
        deltaZ = delta1Z/6 + delta2Z/3 + delta3Z/3 + delta4Z/6
        x += h
        y += deltaY
        z += deltaZ
    return (x,y,z)

lb = 0
ub = duracaoTratamento
pi = (0, 0, 0)
h = 0.1


"""
x = np.arange(0, freq * 5, 1)     # range(0, 90, 1)
y = [ D(i) for i in x ]   # list comprehension
plt.plot(x,y)
plt.show()"""

# MASSA DE AMOXICILINA NO PLASMA

x = np.arange(0, duracaoTratamento, h)     # range(0, 90, 1)
#y = [ rk4(lb, i , h, pi, dmidt, dmpdt)[2] for i in x ]   # list comprehension     
miArray = []                
mpArray = []                                
rk4(lb, duracaoTratamento, h, pi, dmidt, dmpdt)        # Much faster to call the function once and append all y to an array

# plotting the points  
plt.plot(x,mpArray)

# naming the x axis 
plt.xlabel('time (min)') 
# naming the y axis 
plt.ylabel('massa de amoxicilina no plasma (mg) ') 
  
# giving a title to my graph 
plt.title('test') 
  
# function to show the plot 
plt.show()


