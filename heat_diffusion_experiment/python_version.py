from six.moves import xrange
import numpy as np

def python_version(M, total_time=500.0):
    
    # Physical properties --------------------------------------------------------------------------------
    k = 385.0
    rho = 8000.0
    cp = 400.0;
    
    # Geometrical properties -----------------------------------------------------------------------------
    L = 1.0
    l = L/M
    A = l**2
    V = l**3.0        
    
    # Boundary conditions --------------------------------------------------------------------------------
    u = np.zeros((M,M),dtype=np.float64) # Initial temperature
    u[M-1,:] = 25 # Temperature at top wall
    
    # , 'necessary dt: %f'%((0.9 * rho * cp * (l ** 2) / (4.0 * k)))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))

    f = k * A / l # Diffusion term
    B = rho * cp * V / dt # Transient term

    num_timesteps = int(total_time/dt)
    for t in xrange(num_timesteps):
        u0 = np.copy(u)        
        for i in xrange(1,M-1):
            for j in xrange(1,M-1):
                u[i,j] = (f*u0[i,j+1] + f*u0[i,j-1] + f*u0[i+1,j] + f*u0[i-1,j] + B*u0[i,j]) / (B + 4*f)
    return u
