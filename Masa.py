import numpy as np
import matplotlib.pyplot as plt

M0 = 128.5 #fg (masa inicial)
M02 = 128.25 #fg (masa luego de la desorción de moléculas)
DM = 0.23  #fg (masa de 1ML)

m = 0.01 #masa de las moléculas
Rad = 0.043/m

#ecuación para la adsorción
def dMdt(t,M):
    return m*Rad*(1-(M-M02)/DM)

#ecuación para la desorción
def dMdt2(t,M2):
    return -m*Rad*(1+(M2-M0)/DM)

N = 3000
t = np.linspace(1,10,N)     #tiempo donde ocurre la adsorción (láser apagado)
t2 = np.linspace(0,0.8,N)   #tiempo donde ocurre la desorción (láser encendido)
h = 1/100                   #paso
M = np.zeros(t.size, dtype = 'complex_')    #masa durante ads
M2 = np.zeros(t2.size, dtype = 'complex_')  #masa durante des


M2[0] = M0      #masa inicial

#Runge-kutta orden 4
for i in range(t2.size-1):
    k1 = dMdt2(t2[i],M2[i])
    k2 = dMdt2(t2[i]+h/2, M2[i]+h*k1/2)
    k3 = dMdt2(t2[i]+h/2, M2[i]+h*k2/2)
    k4 = dMdt2(t2[i], M2[i]+h*k3)    
    
    M2[i+1] = M2[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

M[0] = M2[N-1] #la adsorción comienza en el punto donde terminó la desorción

for i in range(t.size-1):
    k1 = dMdt(t[i],M[i])
    k2 = dMdt(t[i]+h/2, M[i]+h*k1/2)
    k3 = dMdt(t[i]+h/2, M[i]+h*k2/2)
    k4 = dMdt(t[i], M[i]+h*k3)
    
    M[i+1] = M[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

#plot
plt.figure(figsize=(6.5,5))
plt.rcParams.update({'font.size': 15})
plt.plot(t2, M2-M0,label="Láser encendido") 
plt.plot(t, M-M0,label="Láser apagado")
plt.legend()
plt.xlabel("Horas")
plt.ylabel("Masa-128.5 [fg]")
#plt.savefig("thermal-treatment.png")
