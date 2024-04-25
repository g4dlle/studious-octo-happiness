from scipy.interpolate import interp1d
from constants import *
global Di, R1, R2, R3, R4, R5

Te = 1.2
Di = np.zeros((N,2))
R1 = np.zeros((N,2))
R2 = 6.2e-10
R3 = 2e-7*np.exp(-6.2/(Te))
R4 = 2.7e-13*(Te)**(-3/4)
R5 = 8.75e-39*(Te)**(-9/2)
R6 = 10e-9*(Te)*0.5*np.exp(-11.6/Te)
R7 = 2.5e-11
R8 = 3e-15
R9 = 2e-7
patternD = 'Diffusion coefficient'
patternR = 'C1    Ar    Effective (momentum)'
with open('output.dat', 'r') as file:
    for index, line in enumerate(file):
        if patternD in line:
            indexD = index
        elif patternR in line:
            indexR = index
            break
file = open('output.dat', 'r')
lines = file.readlines()
for i in range(N):
    Di[i] = [float(lines[indexD+1+i][:5]),float(lines[indexD+1+i][-11:-1])]
    R1[i] = [float(lines[indexR+2+i][:5]),float(lines[indexR+2+i][-11:-1])]
file.close()

R1_linear = interp1d(list(i[0] for i in R1), list(i[1] for i in R1), kind = 'linear')
Di_linear = interp1d(list(i[0] for i in Di), list(i[1] for i in Di), kind = 'linear')
