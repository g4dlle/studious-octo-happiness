import numpy as np
from constants import *
from rate_coefficient import *

def EEEE(phi):
    r = np.linspace(hr/2,R,N)
    n1 = len(r)
    n2 = len(z)
    h1 = r[1]-r[0]
    h2 = z[1]-z[0]
    E_r = np.ones((n1, n2), 'float')*10.
    E_z = np.ones((n1, n2), 'float')*10.
    for j in range(0,n2-1):
        E_r[0,j] = (phi[1,j]-phi[0,j])/h1
        for i in range(0,n1-1):
            E_r[i,j] = (phi[i,j]-phi[i+1,j])/h1
    for i in range(0,n1-1):
        E_z[i,0] = (phi[i,1]-phi[i,0])/h2
        for j in range(0,n2-1):
            E_z[i,j] = (phi[i,j]-phi[i,j+1])/h2
    return abs(np.sqrt(E_z**2+E_r**2))

def Rg_pr(FF):
    f = np.zeros((N,M))
    variable_phi=EEEE(FF)/Net
    for i in range(N):
        for j in range(M):
            if variable_phi[i,j] < 1: variable_phi[i,j] = 1
            f[i,j] = R1_linear(float(variable_phi[i,j]))*ne[i,j]*Net[i,j] + R2*nm[i,j]**2 + R3*nm[i,j]*ne[i,j] - R4*ne[i,j]*ni[i,j] - R5*ne[i,j]**2*ni[i,j]
    return f

def Rc_di(FF):
    f = np.zeros((N,M))
    variable_phi=EEEE(FF)/Net
    for i in range(N):
        for j in range(M):
            if variable_phi[i,j] < 1: variable_phi[i,j] = 1
            f[i,j] = Di_linear(float(variable_phi[i,j]))
    return f

def Meta():
    Meta = np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            Meta[i,j] = R6*Net[i,j]*ne[i,j] - R2*nm[i,j]**2 - R3*nm[i,j]*ne[i,j] - R7*nm[i,j] - R8**Net[i,j]*nm[i,j] - R9*nm[i,j]*ne[i,j]
    return Meta