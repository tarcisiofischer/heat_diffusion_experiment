import numpy as np



def solve(
    n_elements,
    total_simulation_time,
    # Physical properties -------------------------------------------------------------------------
    k,
    rho,
    cp,

    # Boundary conditions -------------------------------------------------------------------------
    top_temperature,
    bottom_temperature,
    left_temperature,
    right_temperature,

    after_timestep_callback=None,
):
    '''
    Solve the 2d heat transfer using FEM (explicit).

    :param int n_elements:
        Number of elements in each axis.

    :param function after_timestep_callback:
        A function to be ran after each complete timestep, with signature:
        function(solution_vector, current_timestep)
    '''

    # Geometric properties ------------------------------------------------------------------------
    L = 1.0
    l = L / n_elements
    A = l ** 2
    V = l ** 3.0        

    # Boundary conditions -------------------------------------------------------------------------
    u = np.zeros((n_elements, n_elements), dtype=np.float64) # Initial temperature
    u[n_elements - 1, :] = bottom_temperature
    u[0, :] = top_temperature
    u[:, n_elements - 1] = left_temperature
    u[:, 0] = right_temperature

    # dt condition: dt <= (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    num_timesteps = int(total_simulation_time / dt)

    f = k * A / l # Diffusion term
    B = rho * cp * V / dt # Transient term

    for t in range(num_timesteps):
        u0 = np.copy(u)        
        for i in range(1, n_elements - 1):
            for j in range(1, n_elements - 1):
                u[i, j] = (
                    f * u0[i, j + 1] +
                    f * u0[i, j - 1] +
                    f * u0[i + 1, j] +
                    f * u0[i - 1, j] +
                    B * u0[i, j]
                ) / (B + 4 * f)
        if after_timestep_callback:
            after_timestep_callback(u, t * dt)

    return u
