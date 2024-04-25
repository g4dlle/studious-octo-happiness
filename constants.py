import numpy as np
global eps0, e, p, NE0, E0, F0, U1, U2, mue, mui, alpha, ne_0, k, T, N, M, r, hr, hz
eps0 = (8.85e-12)
e = (1.6e-19)
p = 1e5
NE0 = 1e18
E0 = 5e3
F0 = 1500
U1 = -100
U2 = 0
mue = 1
mui = 1
alpha = 1
ne_0 = 1e18
k = 1.380649e-23
T = 300

N = 4
M = N

ne = np.zeros((N,M))
De = np.zeros((N,M))

L = 1
R = 1
r1 = R
hr = (1/(N-0.5))
hz = (L/(M-1))
tEnd = 1
tau = 1
 
r = np.linspace(hr/2,R,N)

z = np.linspace(0,L,M)

FF = np.zeros((N,M))
Net = np.zeros((N,M))
for j in range(0,M):
    for i in range(0,N):
        FF[i,j] = (-50*z[j])
        Net[i,j] = p/(k*T)
        ne[i][j] = ne_0*(1-(r[i]/R)**2)*(z[j]*(L-z[j]))

ni = ne.copy()
nm = ne.copy()

for i in range(N):
    FF[i,0] = 0 #низ
for i in range(N):
    FF[i,M-1] = -50