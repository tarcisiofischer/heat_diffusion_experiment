import networkx as nx


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
    G.data['k'] = physical_properties.k
    G.data['rho'] = physical_properties.rho
    G.data['c_p'] = physical_properties.c_p

    for node_id in G.nodes():
        node = G.node[node_id]

        # Problem variables
        node['T'] = initial_condition.T

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

    current_time = 0.0
    # solution = initial_condition
    while current_time < timestep_properties.final_time:
        # Aliases
        dx = graph.data['dx']
        dy = graph.data['dy']
        k = graph.data['k']
        rho = graph.data['rho']
        c_p = graph.data['c_p']
        dt = timestep_properties.delta_t

        # Build the coeficient matrix A
        A_P = rho * dx * dy / dt + \
            (2 * k * dy) / (c_p * dx) + \
            (2 * k * dx) / (c_p * dy)
        A_E = (k * dy) / (c_p * dx)
        A_W = (k * dy) / (c_p * dx)
        A_N = (k * dx) / (c_p * dy)
        A_S = (k * dx) / (c_p * dy)

        # Build the matrix B
        # TODO: T_P_old must be get from each node. B_T may be inside the loop or done with
        # numpy array (Probably the first - I wouldn't bother with numpy for now).
        B_T = (rho * dx * dy * T_P_old) / dt
        # for node in G:
        #   build equation for node[i]
        #   append equation to the matrixes

        # Solve Ax = B
        # solve linear system

        # Retrieve the solution
        # update graph

        # Advance in time
        current_time += timestep_properties.delta_t

    return graph
