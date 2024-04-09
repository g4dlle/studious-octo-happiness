# %%
import numpy as np

from parabolic2D import parabolic2D
from constants import *
from rate_coefficient import *
from fsolve import fsolve

# %%
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
    return 

# %%
def Rg_pr(FF):
    f = np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            variable_phi=FF[i,j]/Net[i,j]
            if variable_phi < 1: variable_phi = 1
            f[i,j] = R1_linear(float(variable_phi))*ne[i,j]*Net[i,j] + R2*nm[i,j]**2 + R3*nm[i,j]*ne[i,j] - R4*ne[i,j]*ni[i,j] - R5*ne[i,j]**2*ni[i,j]

# %%
print('решение')
t = 0
while t<= tEnd:
    
    FF,nv = fsolve(FF,ne,ni)
    ne = parabolic2D(r,z,Rg_pr(FF),ne,-1,mue,FF,Di,tau)
    ni = parabolic2D(r,z,Rg_pr(FF),ni,1,mui,FF,Di,tau)
    #nm = parabolic2D(r,z,_,nm,0,mui,FF,Di,tau)
    t = tau+t
    

print((FF.transpose()))
