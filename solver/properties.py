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
    def __init__(self, T, n_x, n_y):
        self.T = np.ones(shape=(n_y, n_x)) * T


class ConstantMapInitialCondition:
    def __init__(self, T):
        self.T = T


class FromFileInitialCondition:
    def __init__(self, filename):
        self.T = np.loadtxt(filename)


class TemperatureBoundaryConditions:
    def __init__(self, east_bc, west_bc, north_bc, south_bc):
        self.east_bc = east_bc
        self.west_bc = west_bc
        self.north_bc = north_bc
        self.south_bc = south_bc


class PrescribedTemperatureBoundaryCondition:
    def __init__(self, T):
        self.T = T


    def type(self):
        from solver.residual_function import BC_TYPE_PRESCRIBED_PHI
        return BC_TYPE_PRESCRIBED_PHI


class PrescribedFlowBoundaryCondition:
    def __init__(self, flow):
        self.prescribed_flow = flow


    def type(self):
        from solver.residual_function import BC_TYPE_PRESCRIBED_FLOW
        return BC_TYPE_PRESCRIBED_FLOW


class TimestepProperties:
    def __init__(self, delta_t, final_time):
        self.delta_t = delta_t
        self.final_time = final_time
