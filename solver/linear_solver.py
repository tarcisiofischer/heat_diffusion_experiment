from copy import deepcopy
import numpy as np
from petsc4py import PETSc



class GeometricProperties:
    def __init__(self, n_x, n_y, size_x, size_y):
        self.n_x = n_x
        self.n_y = n_y
        self.size_x = size_x
        self.size_y = size_y


class PhysicalProperties:
    def __init__(self, k, rho, c_p):
        self.k = k
        self.rho = rho
        self.c_p = c_p


class ConstantInitialCondition:
    def __init__(self, T):
        self.T = T


class GhostNodeBoundaryCondition:
    def __init__(self, T_E, T_W, T_N, T_S):
        self.T_E = T_E
        self.T_W = T_W
        self.T_N = T_N
        self.T_S = T_S


class TimestepProperties:
    def __init__(self, delta_t, final_time):
        self.delta_t = delta_t
        self.final_time = final_time


def build_grid(geometric_properties, physical_properties, initial_condition, boundary_condition):
    # Aliases
    size_x = geometric_properties.size_x
    n_x = geometric_properties.n_x
    size_y = geometric_properties.size_y
    n_y = geometric_properties.n_y

    # Computed variables
    dx = size_x / n_x
    dy = size_y / n_y

    # Data among the grid
    G = {}
    G['dx'] = dx
    G['dy'] = dy
    G['n_x'] = n_x
    G['n_y'] = n_y
    G['k'] = physical_properties.k
    G['rho'] = physical_properties.rho
    G['c_p'] = physical_properties.c_p
    G['T'] = np.ones(shape=(n_x, n_y)) * initial_condition.T

    return G


def transient_term(rho, T_P, T_Po, dx, dy, dt):
    # Accumulative term
    return (rho * T_P - rho * T_Po) * dx * dy / dt


def diffusive_flux_term(k, c_p, T_0, T_1, d0, d1):
    return (k / c_p * (T_0 - T_1) / d0) * d1


def residual_function(snes, X, f, t, graph, old_graph, timestep_properties, boundary_condition):
    n_x = graph['n_x']
    n_y = graph['n_y']
    dt = timestep_properties.delta_t
    dx = graph['dx']
    dy = graph['dy']
    k = graph['k']
    rho = graph['rho']
    c_p = graph['c_p']
    x_old = old_graph['T']

    x = X.getArray(readonly=True).reshape((n_x, n_y))
    eqs = np.zeros(shape=(n_x, n_y))

    # Transient (Accumulative) term
    eqs += transient_term(rho, x[:, :], x_old[:, :], dx, dy, dt)

    # Diffusive term
    # East flux
    eqs[:, :-1] += -diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    eqs[:, -1:] += -diffusive_flux_term(k, c_p, boundary_condition.T_E(t), x[:, -1:], dx / 2.0, dy)
    # West flux
    eqs[:, 1:] += +diffusive_flux_term(k, c_p, x[:, 1:], x[:, :-1], dx, dy)
    eqs[:, :1] += +diffusive_flux_term(k, c_p, x[:, :1], boundary_condition.T_W(t), dx / 2.0, dy)
    # South flux
    eqs[:-1, :] += -diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dx, dy)
    eqs[-1:, :] += -diffusive_flux_term(k, c_p, boundary_condition.T_S(t), x[-1:, :], dx / 2.0, dy)
    # North flux
    eqs[1:, :] += +diffusive_flux_term(k, c_p, x[1:, :], x[:-1, :], dx, dy)
    eqs[:1, :] += +diffusive_flux_term(k, c_p, x[:1, :], boundary_condition.T_N(t), dx / 2.0, dy)

    f[:] = eqs.reshape(n_x * n_y)


def create_solver(use_multigrid, use_newton, n_x, n_y):
    snes = PETSc.SNES().create()

    options = PETSc.Options()
    options.clear()
    if use_multigrid:
        options.setValue('pc_type', 'gamg')
        options.setValue('pc_gamg_reuse_interpolation', 'True')
    else:
        options.setValue('pc_type', 'ilu')
    options.setValue('ksp_type', 'gmres')
    if use_newton:
        options.setValue('snes_type', 'newtonls')
    else:
        options.setValue('snes_type', 'ksponly')
    options.setValue('ksp_atol', 1e-7)
    snes.setFromOptions()

    dmda = PETSc.DMDA().create([n_x, n_y], dof=1, stencil_width=1, stencil_type='star')
    snes.setDM(dmda)
    return snes


def solve(
    geometric_properties,
    physical_properties,
    initial_condition,
    boundary_condition,
    timestep_properties,

    # Output options
    on_timestep_callback=None,

    # Solver options
    use_multigrid=True,
    use_newton=False,
):
    if on_timestep_callback is None:
        on_timestep_callback = lambda *args, **kwargs: None

    graph = build_grid(
        geometric_properties,
        physical_properties,
        initial_condition,
        boundary_condition,
    )

    n_x = geometric_properties.n_x
    n_y = geometric_properties.n_y
    snes = create_solver(use_multigrid, use_newton, n_x, n_y)
    r = PETSc.Vec().createSeq(n_x * n_y)  # residual vector
    x = PETSc.Vec().createSeq(n_x * n_y)  # solution vector
    b = PETSc.Vec().createSeq(n_x * n_y)  # right-hand side

    current_time = 0.0
    old_graph = deepcopy(graph)
    on_timestep_callback(current_time, graph['T'][:])

    # solution = initial_condition
    while current_time < timestep_properties.final_time:
        snes.setFunction(residual_function, r, [current_time, graph, old_graph, timestep_properties, boundary_condition])
        initial_guess = np.array(old_graph['T'])
        x.setArray(initial_guess)
        b.set(0)
        snes.solve(b, x)
        solution = x[:].reshape(n_x, n_y)
        solution.flags['WRITEABLE'] = False

        # Retrieve the solution
        old_graph = deepcopy(graph)
        graph['T'][:] = solution

        # Advance in time
        current_time += timestep_properties.delta_t

        on_timestep_callback(current_time, solution)

    return graph
