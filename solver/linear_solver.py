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


def diffusive_term(k, c_p, T_E, T_P, T_W, T_N, T_S, dx, dy):
    result = 0
    
    # Left side diffusive flux term
    result += +(k / c_p * (T_E - T_P) / dx) * dy
    # Right side diffusive flux term
    result += -(k / c_p * (T_P - T_W) / dx) * dy
    # Top side diffusive flux term
    result += +(k / c_p * (T_N - T_P) / dy) * dx
    # Bottom side diffusive flux term
    result += -(k / c_p * (T_P - T_S) / dy) * dx
    
    return result


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

    if on_timestep_callback is None:
        on_timestep_callback = lambda *args, **kwargs: None

    graph = build_grid(
        geometric_properties,
        physical_properties,
        initial_condition,
        boundary_condition,
    )
    old_graph = deepcopy(graph)

    current_time = 0.0
    on_timestep_callback(current_time, graph['T'][:])

    n_x = graph['n_x']
    n_y = graph['n_y']
    dmda = PETSc.DMDA().create([n_x, n_y], dof=1, stencil_width=1, stencil_type='star')

    dt = timestep_properties.delta_t
    dx = graph['dx']
    dy = graph['dy']
    k = graph['k']
    rho = graph['rho']
    c_p = graph['c_p']

    # Build the matrix coeficients
    A_P = rho * dx * dy / dt + \
        2 * (k * dy) / (c_p * dx) + \
        2 * (k * dx) / (c_p * dy)
    A_E = (k * dy) / (c_p * dx)
    A_W = (k * dy) / (c_p * dx)
    A_N = (k * dx) / (c_p * dy)
    A_S = (k * dx) / (c_p * dy)

    snes = PETSc.SNES().create()
    snes.setFromOptions()

    r = PETSc.Vec().createSeq(n_x * n_y)  # residual vector
    x = PETSc.Vec().createSeq(n_x * n_y)  # solution vector
    b = PETSc.Vec().createSeq(n_x * n_y)  # right-hand side

    # solution = initial_condition
    while current_time < timestep_properties.final_time:
        # Aliases
        T_old = old_graph['T']
        B_T = (rho * dx * dy * T_old) / dt

        # Build the equation system
        def residual_function(snes, X, f):
            x = X.getArray(readonly=True).reshape((n_x, n_y))
            eqs = np.zeros(shape=(n_x, n_y))
            eqs[1:-1, 1:-1] = \
                  transient_term(rho, x[1:-1, 1:-1], T_old[1:-1, 1:-1], dx, dy, dt) \
                - diffusive_term(k, c_p, x[1:-1, 2:], x[1:-1, 1:-1], x[1:-1, :-2], x[:-2, 1:-1], x[2:, 1:-1], dx, dy)

            # Boundary conditions
            # Left
            eqs[1:-1, :1] = \
                  A_P * x[1:-1, :1] \
                - boundary_condition.T_W(current_time) \
                - A_E * x[1:-1, 1:2] \
                - A_N * x[:-2, :1] \
                - A_S * x[2:, :1] \
                - B_T[1:-1, :1]
            # Right
            eqs[1:-1, -1:] = \
                  A_P * x[1:-1, -1:] \
                - A_W * x[1:-1, -2:-1] \
                - boundary_condition.T_E(current_time) \
                - A_N * x[:-2, -1:] \
                - A_S * x[2:, -1:] \
                - B_T[1:-1, -1:]
            # Bottom
            eqs[-1:, 1:-1] = \
                  A_P * x[-1:, 1:-1] \
                - A_W * x[-1:, :-2] \
                - A_E * x[-1:, 2:] \
                - A_N * x[-2:-1, 1:-1] \
                - boundary_condition.T_S(current_time) \
                - B_T[-1:, 1:-1]
            # Top
            eqs[:1, 1:-1] = \
                  A_P * x[:1, 1:-1] \
                - A_W * x[:1, :-2] \
                - A_E * x[:1, 2:] \
                - boundary_condition.T_N(current_time) \
                - A_S * x[1:2, 1:-1] \
                - B_T[:1, 1:-1]

            # Top-Left
            eqs[:1, :1] = \
                  A_P * x[:1, :1] \
                - boundary_condition.T_W(current_time) \
                - A_E * x[:1, 1:2] \
                - boundary_condition.T_N(current_time) \
                - A_S * x[1:2, :1] \
                - B_T[:1, :1]
            # Top-Right
            eqs[:1, -1:] = \
                  A_P * x[:1, -1:] \
                - A_W * x[:1, -2:-1] \
                - boundary_condition.T_E(current_time) \
                - boundary_condition.T_N(current_time) \
                - A_S * x[1:2, -1:] \
                - B_T[:1, -1:]
            # Bottom-Left
            eqs[-1:, :1] = \
                  A_P * x[-1:, :1] \
                - boundary_condition.T_W(current_time) \
                - A_E * x[-1:, 1:2] \
                - A_N * x[-2:-1, :1] \
                - boundary_condition.T_S(current_time) \
                - B_T[-1:, :1]
            # Bottom-Right
            eqs[-1:, -1:] = \
                  A_P * x[-1:, -1:] \
                - A_W * x[-1:, -2:-1] \
                - boundary_condition.T_E(current_time) \
                - A_N * x[-2:-1, -1:] \
                - boundary_condition.T_S(current_time) \
                - B_T[-1:, -1:]

            f[:] = eqs.reshape(n_x * n_y)

        # Solve Ax = B
        snes.setFunction(residual_function, r)

        snes.setDM(dmda)

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
