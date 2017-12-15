import numpy as np



def transient_term(rho, T_P, T_Po, dx, dy, dt):
    # Accumulative term
    return (rho * T_P - rho * T_Po) * dx * dy / dt


def diffusive_flux_term(k, c_p, T_0, T_1, d0, d1):
    return (k / c_p * (T_0 - T_1) / d0) * d1


def residual_function(
    snes,
    X,
    f,
    t,
    grid,
    old_grid,
    timestep_properties,
    boundary_conditions
):
    n_x = grid['n_x']
    n_y = grid['n_y']
    dt = timestep_properties.delta_t
    dx = grid['dx']
    dy = grid['dy']
    k = grid['k']
    rho = grid['rho']
    c_p = grid['c_p']
    x_old = old_grid['T']

    x = X.getArray(readonly=True).reshape((n_x, n_y))
    eqs = np.zeros(shape=(n_x, n_y))

    # Transient (Accumulative) term
    eqs += transient_term(rho, x[:, :], x_old[:, :], dx, dy, dt)

    # Diffusive term
    # East flux
    eqs[:, :-1] += -diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    eqs[:, -1:] += -diffusive_flux_term(k, c_p, boundary_conditions.east_bc.T(t), x[:, -1:], dx / 2.0, dy)
    # West flux
    eqs[:, 1:] += +diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    eqs[:, :1] += +diffusive_flux_term(k, c_p, x[:, :1], boundary_conditions.west_bc.T(t), dx / 2.0, dy)
    # South flux
    eqs[:-1, :] += -diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dx, dy)
    eqs[-1:, :] += -diffusive_flux_term(k, c_p, boundary_conditions.south_bc.T(t), x[-1:, :], dx / 2.0, dy)
    # North flux
    eqs[1:, :] += +diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dx, dy)
    eqs[:1, :] += +diffusive_flux_term(k, c_p, x[:1, :], boundary_conditions.north_bc.T(t), dx / 2.0, dy)

    f[:] = eqs.reshape(n_x * n_y)
