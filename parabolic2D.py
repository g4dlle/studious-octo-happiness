import numpy as np
def decLU3(c, a, b):
    # print("Это LU3")
    # print("c=", c,"a=", a,"b=", b)
    """
    Input of tridiagonal matrix A:
        a[i] = A[i,i],
        b[i] = A[i,i+1],
        c[i] = A[i,i-1]
    Returns thr decomposition LU for tridiagonal matrix:
        d[i] = L[i,i]
        u[i] = U[i,i+1]
        l[i] = L[i,i-1]
    """
    n = len(a)
    d = np.copy(a)
    u = np.copy(b)
    l = np.copy(c)
    for i in range(1, n):
        al = l[i] / d[i-1]
        d[i] = d[i] - al*u[i-1]
        l[i] = al
    return d, u, l
def solveLU3(c, a, b, f):
    """
    Solve the linear system Ax = b with tridiagonal matrix:
        a[i] = A[i,i]
        b[i] = A[i,i+1],
        c[i] = A[i,i-1]
    """
    n = len(a)
    # LU decomposition
    d, u, l = decLU3(c, a, b)
    x = np.copy(f)
    # forward substitution process
    for i in range(1, n):
        x[i] = x[i] - l[i]*x[i-1]
    # back substitution process
    x[n-1] = x[n-1] / d[n-1]
    for i in range(n-2,-1,-1):
        x[i] = (x[i] - u[i]*x[i+1]) / d[i]
    return x
def parabolic2D(r, z, f, y_0, znak, mu, phi, D, tau):
    """
    znak -- знак заряда частицы (+1 для иона, -1 для электрона)
    Numerical Solution of the Dirichlet problem
    for two-dimentional parabolic equation.
    Use additive difference scheme of alternating directions.
    r: 1D NumPy array that contains the x-coordinates of the grid points.
    z: 1D NumPy array that contains the y-coordinates of the grid points.
    f: 2D NumPy array that contains the values of the right-hand side function at each grid point.
    y_0: 2D NumPy array that contains the initial guess for the solution.
    znak: scalar that represents the sign of the charge of the particles (1 for an ion, -1 for an electron).
    mu: scalar that represents the magnetic permeability of the medium.
    phi: 2D NumPy array that contains the values of the electric potential at each grid point.
    D: 2D NumPy array that contains the values of the diffusion coefficient in the x and y directions, respectively.
    tau: scalar that represents the time step.
    """
    # print('PARABOLIC2d')
    # print('Нач. приближение \n', y_0)
    n1 = len(r)
    n2 = len(z)
    h1 = r[1]-r[0]
    h2 = z[1]-z[0]
    
    a1 = np.zeros((n1), 'float')
    b1 = np.ones((n1), 'float')
    c1 = np.zeros((n1), 'float')
    q1 = np.zeros((n1), 'float')
    a2 = np.zeros((n2), 'float')
    b2 = np.ones((n2), 'float')
    c2 = np.zeros((n2), 'float')
    q2 = np.zeros((n2), 'float')
    
    D_r = np.ones((n1, n2), 'float')
    E_r = np.ones((n1, n2), 'float')*10.
    alpha_r = np.ones((n1, n2), 'float')
    I_0_r = np.ones((n1, n2), 'float')
    
    D_z = np.ones((n1, n2), 'float')
    E_z = np.ones((n1, n2), 'float')*10.
    alpha_z = np.ones((n1, n2), 'float')
    I_0_z = np.ones((n1, n2), 'float')
    
    y0 = y_0
    y = np.copy(y_0)
    
    # while t0 < tEnd - 0.001*tau:
    # x1 direction
    # вычисляем коэффициенты в промежуточных узлах
    for j in range(0,n2-1):
        D_r[0,j] = (D[0,j]+D[1,j])/2
        E_r[0,j] = (phi[1,j]-phi[0,j])/h1
        alpha_1 = -znak*mu*E_r[0,j]*h1/D_r[0,j]
        alpha_r[0,j] = alpha_1
        if np.isclose(alpha_1, 0.0, atol=1e-8):
            I_0_r[0,j] = 1
        else:
            I_0_r[0,j] = (np.exp(alpha_r[0,j])-1)/alpha_r[0,j]
        
        for i in range(0,n1-1):
            D_r[i,j] = (D[i,j]+D[i+1,j])/2
            E_r[i,j] = (phi[i,j]-phi[i+1,j])/h1
            alpha_1 = -znak*mu*E_r[i,j]*h1/D_r[i,j]
            alpha_r[i,j] = alpha_1
            if np.isclose(alpha_1, 0, atol=1e-8):
                I_0_r[i,j] = 1
            else:
                I_0_r[i,j] = (np.exp(alpha_r[i,j])-1)/alpha_r[i,j]
            
    for i in range(0,n1-1):
        D_z[i,0] = (D[i,0]+D[i,1])/2
        E_z[i,0] = (phi[i,1]-phi[i,0])/h2
        alpha_1 = -znak*mu*E_z[i,0]*h2/D_z[i,0]
        alpha_z[i,0] = alpha_1
        if np.isclose(alpha_1, 0.0, atol=1e-8):
            I_0_z[i,0] = 1.
        else:
            I_0_z[i,0] = (np.exp(alpha_z[i,0])-1)/alpha_z[i,0]
        
        for j in range(0,n2-1):
            D_z[i,j] = (D[i,j]+D[i,j+1])/2
            E_z[i,j] = (phi[i,j]-phi[i,j+1])/h2
            alpha_1 = -znak*mu*E_z[i,j]*h2/D_z[i,j]
            alpha_z[i,j] = alpha_1
            if np.isclose(alpha_1, 0.0, atol=1e-8):
                I_0_z[i,j] = 1.
            else:
                I_0_z[i,j] = (np.exp(alpha_z[i,j])-1)/alpha_z[i,j]                

    for j in range(1,n2-1):
        b1[0] = (2/h1**2)*D_r[0,j]/I_0_r[0,j] + 2./tau
        c1[0] = -(2/h1**2)*D_r[0,j]/I_0_r[0,j]
        q1[0] = (0.5*f[0,j]) + 2*y[0,j]/tau \
            + (D_z[0,j-1]/(I_0_z[0,j-1])*y[0,j-1] \
                - (D_z[0,j-1]/(I_0_z[0,j-1])*np.exp(alpha_z[0,j-1]) \
                  + D_z[0,j]/(I_0_z[0,j]))*y[0,j] \
                    + D_z[0,j]/(I_0_z[0,j])*np.exp(alpha_z[0,j])*y[0,j+1])/h2**2
        
        for i in range(1,n1-1):
            a1[i] = -((r[i]-h1/2))*D_r[i-1,j]/(r[i]*h1**2*I_0_r[i-1,j])
            b1[i] = ((r[i]-h1/2))*D_r[i-1,j]/(r[i]*h1**2*I_0_r[i-1,j])\
                *np.exp(alpha_r[i-1,j]) \
                + (r[i]+h1/2)/r[i]*D_r[i,j]/(h1**2*I_0_r[i,j]) + 2. / tau
            c1[i] = -(r[i]+h1/2)/r[i]*D_r[i,j]/(h1**2*I_0_r[i,j]) \
                *np.exp(alpha_r[i,j])
            q1[i] = f[i,j]/2. + 2.*y[i,j]/tau \
                + (D_z[i,j-1]/(h2*I_0_z[i,j-1])*y[i,j-1] \
                   - (D_z[i,j-1]/(h2*I_0_z[i,j-1])*np.exp(alpha_z[i,j-1]) \
                      + D_z[i,j]/(h2*I_0_z[i,j]))*y[i,j] \
                       + D_z[i,j]/(h2*I_0_z[i,j])*np.exp(alpha_z[i,j])*y[i,j+1])/h2   
        y0[:,j] = solveLU3(a1, b1, c1, q1)
        
    for j in range(1,n2-1):
        a2[j] = -D_z[0,j-1]/(h2**2*I_0_z[0,j-1])
        b2[j] = -a2[j]*np.exp(alpha_z[0,j-1])+D_z[0,j]/(h2**2*I_0_z[0,j]) \
            + 2./tau
        c2[j] = -D_z[0,j]/(h2**2*I_0_z[0,j])*np.exp(alpha_z[0,j])
        q2[j] = f[0, j]/2. + 2.*y0[0,j]/tau \
            - (((r[0]+h1/2)/r[0]*D_r[0,j]/(I_0_r[0,j]))*y0[0,j]\
                   - (r[0]+h1/2)/r[0]*D_r[0,j]/(I_0_r[0,j])\
                       *np.exp(alpha_r[0,j])*y0[1,j])/h1**2    
    y[0,:] = solveLU3(a2, b2, c2, q2)
    
    # x2 direction
    for i in range(1,n1-1):
        for j in range(1,n2-1):
            a2[j] = -D_z[i,j-1]/(h2**2*I_0_z[i,j-1])
            b2[j] = -a2[j]*np.exp(alpha_z[i,j-1])+D_z[i,j]/(h2**2*I_0_z[i,j]) \
                + 2. / tau
            c2[j] = -D_z[i,j]/(h2**2*I_0_z[i,j])*np.exp(alpha_z[i,j])
            q2[j] = f[i, j]/2. + 2.*y0[i,j]/tau\
                - (-(r[i]-h1/2)/D_r[i-1,j]/(I_0_r[i-1,j])*y0[i-1,j]\
                   + ((r[i]-h1/2)/D_r[i-1,j]/(I_0_r[i-1,j])\
                      *np.exp(alpha_r[i-1,j])\
                      + (r[i]+h1/2)/D_r[i,j]/(I_0_r[i,j]))*y0[i,j]\
                       - (r[i]+h1/2)/D_r[i,j]/(I_0_r[i,j])\
                           *np.exp(alpha_r[i,j])*y0[i+1,j])/r[i]/h1**2      
        y[i,:] = solveLU3(a2, b2, c2, q2)
    return y