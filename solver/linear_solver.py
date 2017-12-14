import networkx as nx
from copy import deepcopy
import numpy as np



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


def build_graph(geometric_properties, physical_properties, initial_condition, boundary_condition):
    '''
    Builds the problem's non-oriented graph G, where each node represents a control volume, and each
    edge represent that two control volumes are neighbors.
    '''
    # Aliases
    size_x = geometric_properties.size_x
    n_x = geometric_properties.n_x
    size_y = geometric_properties.size_y
    n_y = geometric_properties.n_y

    # Computed variables
    dx = size_x / n_x
    dy = size_y / n_y

    G = nx.grid_2d_graph(n_x, n_y)

    # Constants data among the grid
    G.data = {}
    G.data['dx'] = dx
    G.data['dy'] = dy
    G.data['n_x'] = n_x
    G.data['n_y'] = n_y
    G.data['k'] = physical_properties.k
    G.data['rho'] = physical_properties.rho
    G.data['c_p'] = physical_properties.c_p
    G.data['T'] = np.ones(shape=(n_x, n_y)) * initial_condition.T

    return G


def build_matrix_structure(graph):
    return None


def solve(
    geometric_properties,
    physical_properties,
    initial_condition,
    boundary_condition,
    timestep_properties,
):
    graph = build_graph(
        geometric_properties,
        physical_properties,
        initial_condition,
        boundary_condition,
    )
    old_graph = deepcopy(graph)

    current_time = 0.0
    # solution = initial_condition
    while current_time < timestep_properties.final_time:
        # Aliases
        dt = timestep_properties.delta_t
        dx = graph.data['dx']
        dy = graph.data['dy']
        k = graph.data['k']
        rho = graph.data['rho']
        c_p = graph.data['c_p']
        n_x = graph.data['n_x']
        n_y = graph.data['n_y']

        # Build the matrix coeficients
        A_P = rho * dx * dy / dt + \
            2 * (k * dy) / (c_p * dx) + \
            2 * (k * dx) / (c_p * dy)
        A_E = (k * dy) / (c_p * dx)
        A_W = (k * dy) / (c_p * dx)
        A_N = (k * dx) / (c_p * dy)
        A_S = (k * dx) / (c_p * dy)
        B_T = (rho * dx * dy * old_graph.data['T']) / dt

        # Build the equation system
        def residual_function(snes, X, f):
            x = np.array(X).reshape((n_x, n_y))
            eqs = np.zeros(shape=(n_x, n_y))
            eqs[1:-1, 1:-1] = \
                  A_P * x[1:-1, 1:-1] \
                - A_W * x[1:-1, :-2] \
                - A_E * x[1:-1, 2:] \
                - A_N * x[:-2, 1:-1] \
                - A_S * x[2:, 1:-1] \
                - B_T[1:-1, 1:-1]

            # Boundary conditions
            # Left
            eqs[1:-1, :1] = \
                  A_P * x[1:-1, :1] \
                - A_W * boundary_condition.T_W \
                - A_E * x[1:-1, 1:2] \
                - A_N * x[:-2, :1] \
                - A_S * x[2:, :1] \
                - B_T[1:-1, :1]
            # Right
            eqs[1:-1, -1:] = \
                  A_P * x[1:-1, -1:] \
                - A_W * x[1:-1, -2:-1] \
                - A_E * boundary_condition.T_E \
                - A_N * x[:-2, -1:] \
                - A_S * x[2:, -1:] \
                - B_T[1:-1, -1:]
            # Bottom
            eqs[-1:, 1:-1] = \
                  A_P * x[-1:, 1:-1] \
                - A_W * x[-1:, :-2] \
                - A_E * x[-1:, 2:] \
                - A_N * x[-2:-1, 1:-1] \
                - A_S * boundary_condition.T_S \
                - B_T[-1:, 1:-1]
            # Top
            eqs[:1, 1:-1] = \
                  A_P * x[:1, 1:-1] \
                - A_W * x[:1, :-2] \
                - A_E * x[:1, 2:] \
                - A_N * boundary_condition.T_N \
                - A_S * x[1:2, 1:-1] \
                - B_T[:1, 1:-1]

            # Top-Left
            eqs[:1, :1] = \
                  A_P * x[:1, :1] \
                - A_W * boundary_condition.T_W \
                - A_E * x[:1, 1:2] \
                - A_N * boundary_condition.T_N \
                - A_S * x[1:2, :1] \
                - B_T[:1, :1]
            # Top-Right
            eqs[:1, -1:] = \
                  A_P * x[:1, -1:] \
                - A_W * x[:1, -2:-1] \
                - A_E * boundary_condition.T_E \
                - A_N * boundary_condition.T_N \
                - A_S * x[1:2, -1:] \
                - B_T[:1, -1:]
            # Bottom-Left
            eqs[-1:, :1] = \
                  A_P * x[-1:, :1] \
                - A_W * boundary_condition.T_W \
                - A_E * x[-1:, 1:2] \
                - A_N * x[-2:-1, :1] \
                - A_S * boundary_condition.T_S \
                - B_T[-1:, :1]
            # Bottom-Right
            eqs[-1:, -1:] = \
                  A_P * x[-1:, -1:] \
                - A_W * x[-1:, -2:-1] \
                - A_E * boundary_condition.T_E \
                - A_N * x[-2:-1, -1:] \
                - A_S * boundary_condition.T_S \
                - B_T[-1:, -1:]

            f[:] = eqs.reshape(n_x * n_y)

        # Solve Ax = B
        from petsc4py import PETSc
        snes = PETSc.SNES().create()
        r = PETSc.Vec().createSeq(n_x * n_y)  # residual vector
        x = PETSc.Vec().createSeq(n_x * n_y)  # solution vector
        b = PETSc.Vec().createSeq(n_x * n_y)  # right-hand side
        snes.setFunction(residual_function, r)

#         dmda = PETSc.DMDA().create([n_x * n_y], dof=1, stencil_width=n_x*n_y, stencil_type='star')
#         snes.setDM(dmda)

        from scipy.sparse.csr import csr_matrix
        dm = PETSc.DMShell().create()
        j_structure = csr_matrix(np.ones(shape=(n_x*n_y, n_x*n_y)))
        csr = (j_structure.indptr, j_structure.indices, j_structure.data)
        jac = PETSc.Mat().createAIJWithArrays(j_structure.shape, csr)
        dm.setMatrix(jac)
        snes.setDM(dm)

        initial_guess = np.array(old_graph.data['T'])
        x.setArray(initial_guess)
        b.set(0)
        snes.solve(b, x)
        solution = np.array(x).reshape(n_x, n_y)

        # Retrieve the solution
        old_graph = deepcopy(graph)
        graph.data['T'][:] = solution

        # Advance in time
        current_time += timestep_properties.delta_t
        print(current_time)

    import matplotlib.pyplot as plt
    plt.imshow(graph.data['T'])
    plt.show()

    return graph
