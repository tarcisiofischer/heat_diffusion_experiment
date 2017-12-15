import numpy as np



BC_TYPE_PRESCRIBED_PHI = 'prescribed_phi'
BC_TYPE_PRESCRIBED_FLOW = 'prescribed_flow'


def accumulation_term(rho, T_P, T_Po, dx, dy, dt):
    return (rho * T_P - rho * T_Po) * dx * dy / dt


def diffusive_flux_term(k, c_p, T_0, T_1, d0, d1):
    return (k / c_p * (T_0 - T_1) / d0) * d1


def add_transient_term(eqs, rho, x, x_old, dx, dy, dt):
    eqs += accumulation_term(rho, x, x_old, dx, dy, dt)


def add_diffusive_term(eqs, boundary_conditions, n_x, n_y, k, c_p, x, dx, dy, t):
    # East flux
    eqs[:, :-1] += -diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    bc = boundary_conditions.east_bc
    if bc.type() == BC_TYPE_PRESCRIBED_PHI:
        eqs[:, -1:] += -diffusive_flux_term(k, c_p, bc.T(t), x[:, -1:], dx / 2.0, dy)
    elif bc.type() == BC_TYPE_PRESCRIBED_FLOW:
        eqs[:, -1:] += -bc.prescribed_flow(t)
    else:
        assert False, "Unknown boundary condition of type %s" % (bc.type(),)
    # West flux
    eqs[:, 1:] += +diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    bc = boundary_conditions.west_bc
    if bc.type() == BC_TYPE_PRESCRIBED_PHI:
        eqs[:, :1] += +diffusive_flux_term(k, c_p, x[:, :1], bc.T(t), dx / 2.0, dy)
    elif bc.type() == BC_TYPE_PRESCRIBED_FLOW:
        eqs[:, :1] += -bc.prescribed_flow(t)
    else:
        assert False, "Unknown boundary condition of type %s" % (bc.type(),)
    # South flux
    eqs[:-1, :] += -diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dy, dx)
    bc = boundary_conditions.south_bc
    if bc.type() == BC_TYPE_PRESCRIBED_PHI:
        eqs[-1:, :] += -diffusive_flux_term( k, c_p, bc.T(t), x[-1:, :], dy / 2.0, dx)
    elif bc.type() == BC_TYPE_PRESCRIBED_FLOW:
        eqs[-1:, :] += -bc.prescribed_flow(t)
    else:
        assert False, "Unknown boundary condition of type %s" % (bc.type(),)
    # North flux
    eqs[1:, :] += +diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dy, dx)
    bc = boundary_conditions.north_bc
    if bc.type() == BC_TYPE_PRESCRIBED_PHI:
        eqs[:1, :] += +diffusive_flux_term(k, c_p, x[:1, :], bc.T(t), dy / 2.0, dx)
    elif bc.type() == BC_TYPE_PRESCRIBED_FLOW:
        eqs[:1, :] += -bc.prescribed_flow(t)
    else:
        assert False, "Unknown boundary condition of type %s" % (bc.type(),)


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

    add_transient_term(eqs, rho, x, x_old, dx, dy, dt)
    add_diffusive_term(eqs, boundary_conditions, n_x, n_y, k, c_p, x, dx, dy, t)

    f[:] = eqs.reshape(n_x * n_y)
