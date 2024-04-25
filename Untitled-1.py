import numpy as np
from parabolic2D import parabolic2D
from constants import *
from rate_coefficient import *
from fsolve import fsolve
from Right_side import Rg_pr, Rc_di, Meta

print('решение')
t = 0
while t<= tEnd:
    
    FF,nv = fsolve(FF,ne,ni)
    ne = parabolic2D(r,z,Rg_pr(FF),ne,-1,mue,FF,Rc_di(FF),tau)
    ni = parabolic2D(r,z,Rg_pr(FF),ni,1,mui,FF,Rc_di(FF),tau)
    nm = parabolic2D(r,z,Meta(),nm,0,mui,FF,Rc_di(FF),tau)
    t = tau+t

print((FF.transpose()))