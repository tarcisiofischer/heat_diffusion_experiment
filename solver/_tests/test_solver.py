from solver.linear_solver import solve, GeometricProperties, PhysicalProperties, \
    ConstantInitialCondition, TimestepProperties, GhostNodeBoundaryCondition


def test_solver():
    solve(
        GeometricProperties(
            n_x=10,
            n_y=10,
            size_x=1.0,
            size_y=1.0,
        ),
        PhysicalProperties(
            k=1.0,
            rho=1.0,
            c_p=1.0,
        ),
        ConstantInitialCondition(
            T=10.0,
        ),
        GhostNodeBoundaryCondition(
            T_E=11.0,
            T_W=15.0,
            T_N=20.0,
            T_S=30.0,
        ),
        TimestepProperties(
            delta_t=0.1,
            final_time=10.0,
        )
    )
