from solver.linear_solver import GeometricProperties, build_graph, PhysicalProperties, \
    ConstantInitialCondition, GhostNodeBoundaryCondition
import matplotlib.pyplot as plt
import networkx as nx


def test_build_graph():
    g = build_graph(
        GeometricProperties(
            n_x=10,
            n_y=10,
            size_x=100,
            size_y=100,
        ),
        PhysicalProperties(
            k=1.0,
            rho=1.0,
            cp=1.0,
        ),
        ConstantInitialCondition(
            T=10.0,
        ),
        GhostNodeBoundaryCondition(
            T_E=10.0,
            T_W=10.0,
            T_N=10.0,
            T_S=10.0,
        ),
    )
