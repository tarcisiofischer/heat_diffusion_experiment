import numpy as np



def solve(n_elements, num_timesteps=0, total_simulation_time=500.0):
    '''
    Solve the 2d heat transfer using FEM (explicit).

    :param int n_elements:
        Number of elements in each axis.

    :param int num_timesteps:
        If None, will try to calculate based on dt.
    
    :param float total_simulation_time:
        The total simulation time in seconds.
    '''
    # Physical properties -------------------------------------------------------------------------
    k = 385.0
    rho = 8000.0
    cp = 400.0

    # Geometric properties ------------------------------------------------------------------------
    L = 1.0
    l = L / n_elements
    A = l ** 2
    V = l ** 3.0        

    # Boundary conditions -------------------------------------------------------------------------
    u = np.zeros((n_elements, n_elements), dtype=np.float64) # Initial temperature
    u[n_elements - 1, :] = 25 # Temperature at top wall

    # dt condition: dt <= (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    if num_timesteps == 0:
        num_timesteps = int(total_simulation_time / dt)
    
    f = k * A / l # Diffusion term
    B = rho * cp * V / dt # Transient term

    for t in xrange(num_timesteps):
        u0 = np.copy(u)        
        for i in xrange(1, n_elements - 1):
            for j in xrange(1, n_elements - 1):
                u[i, j] = (
                    f * u0[i, j + 1] +
                    f * u0[i, j - 1] +
                    f * u0[i + 1, j] +
                    f * u0[i - 1, j] +
                    B * u0[i, j]
                ) / (B + 4 * f)
    return u
