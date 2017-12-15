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


class PrescribedTemperatureBoundaryCondition:
    def __init__(self, T_E, T_W, T_N, T_S):
        self.T_E = T_E
        self.T_W = T_W
        self.T_N = T_N
        self.T_S = T_S


class TimestepProperties:
    def __init__(self, delta_t, final_time):
        self.delta_t = delta_t
        self.final_time = final_time
