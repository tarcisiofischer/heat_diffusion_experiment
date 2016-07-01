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
    :param int n_elements:
        Number of elements in each axis.
    
    :param float total_simulation_time:
        Total amount of simulation time in seconds
    
    :param float k:
    :param float rho:
    :param float cp:
    :param float top_temperature:
    :param float bottom_temperature:
    :param float left_temperature:
    :param float right_temperature:

    :param function after_timestep_callback:
        A function to be ran after each complete timestep, with signature:
        function(solution_vector, current_timestep)
    '''

    # Geometric properties ------------------------------------------------------------------------
    L = 1.0
    l = L / n_elements

    # Initial condition ---------------------------------------------------------------------------
    u = np.zeros((n_elements, n_elements), dtype=np.float64) # Initial temperature

    # dt condition: dt <= (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    dt = (0.9 * rho * cp * (l ** 2) / (4.0 * k))
    num_timesteps = int(total_simulation_time / dt)

    f = k / l # Diffusion term
    B = rho * cp * l / dt # Transient term

    for t in range(num_timesteps):
        u0 = np.copy(u)

        for i in range(0, n_elements):
            for j in range(0, n_elements):
                u[i, j] = 0.0

                if j + 1 <= n_elements - 1:
                    u[i, j] += f * u0[i, j + 1]
                else:
                    u[i, j] += f / 2.0 * right_temperature

                if j - 1 >= 0:
                    u[i, j] += f * u0[i, j - 1]
                else:
                    u[i, j] += f / 2.0 * left_temperature

                if i + 1 <= n_elements - 1:
                    u[i, j] += f * u0[i + 1, j]
                else:
                    u[i, j] += f / 2.0 * top_temperature

                if i - 1 >= 0:
                    u[i, j] += f * u0[i - 1, j]
                else:
                    u[i, j] += f / 2.0 * bottom_temperature

                u[i, j] += B * u0[i, j]
                u[i, j] /= B + 4 * f
        if after_timestep_callback:
            after_timestep_callback(u, t * dt)

    return u
