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


class TemperatureBoundaryConditions:
    def __init__(self, east_bc, west_bc, north_bc, south_bc):
        self.east_bc = east_bc
        self.west_bc = west_bc
        self.north_bc = north_bc
        self.south_bc = south_bc


class PrescribedTemperatureBoundaryCondition:
    def __init__(self, T):
        self.T = T


class TimestepProperties:
    def __init__(self, delta_t, final_time):
        self.delta_t = delta_t
        self.final_time = final_time
