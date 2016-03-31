import numpy as np
cimport numpy as np
ctypedef np.float64_t DTYPE_t
import cython
from cython.parallel import prange, parallel

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.float64_t, ndim=2] cython_version(
    int M, 
    int num_timesteps=0, 
    double total_time=500.0
    ):
    
    cdef double rho
    cdef double L, k, cp, l, A, V, f, B, dt
    cdef int t, i, j
    cdef double[:,::1] u, u0
    
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
    u0 = np.zeros((M,M),dtype=np.float64) # Initial temperature
    for i in range(M):
        u[M-1,i] = 25.0
    
    # , 'necessary dt: %f'%((0.9 * rho * cp * (l ** 2) / (4.0 * k)))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    if num_timesteps == 0:
        num_timesteps = int(total_time/dt)

    f = k * A / l # Diffusion term
    B = rho * cp * V / dt # Transient term
    
    for t in range(num_timesteps):           
        u0 = np.copy(u)        
        for i in range(1,M-1):
            for j in range(1,M-1):
                u[i,j] = (f*u0[i,j+1] + f*u0[i,j-1] + f*u0[i+1,j] + f*u0[i-1,j] + B*u0[i,j]) / (B + 4*f)
    return np.asarray(u) 
 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.float64_t, ndim=2] cython_parallel_version(
    int M, 
    int num_timesteps=0, 
    double total_time=500.0
    ):
     
    cdef double rho
    cdef double L, k, cp, l, A, V, f, B, dt
    cdef int t, i, j
    cdef double[:,::1] u, u0
     
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
    u0[:] = u 
    for i in range(M):
        u[M-1,i] = 25.
     
    # , 'necessary dt: %f'%((0.9 * rho * cp * (l ** 2) / (4.0 * k)))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    if num_timesteps == 0:
        num_timesteps = int(total_time/dt)
    f = k * A / l # Diffusion term
    B = rho * cp * V / dt # Transient term
    
 
    for t in range(num_timesteps):           
        u0[:] = u        
        for i in prange(1, M-1, nogil=True):
            for j in range(1,M-1):
                u[i,j] = (f*u0[i,j+1] + f*u0[i,j-1] + f*u0[i+1,j] + f*u0[i-1,j] + B*u0[i,j]) / (B + 4*f)
    return np.asarray(u) 

def show(a):
    import matplotlib
    import matplotlib.pyplot as plt
    plt.pcolormesh(np.linspace(0, 1, a.shape[0]) , np.linspace(0, 1, a.shape[1]) , a) 
    plt.show()