import numpy as np
from solver.nonlinear_solver_wrappers import scipy_wrapper,\
    petsc4py_wrapper
from solver.time_stepper import transient_solve
from solver.properties import GeometricProperties, PhysicalProperties,\
    ConstantMapInitialCondition, TemperatureBoundaryConditions,\
    PrescribedFlowBoundaryCondition, TimestepProperties
import pytest

@pytest.mark.parametrize(
    'solver_type',
    (
        'scipy',
        'petsc4py'
    )
)
def test_square_diffusion_smoke(solver_type):
    '''
    Smoke test just to be sure nothing crashes
    '''
    n_x = 100
    n_y = 100
    
    initial_condition = np.zeros(shape=(n_y, n_x))
    initial_condition[
        int(n_y / 2 - .05 * n_y):int(n_y / 2 + .05 * n_y),
        int(n_x / 2 - .05 * n_x):int(n_x / 2 + .05 * n_x)
    ] = 100.0
    
    if solver_type == 'scipy':
        nonlin_package = scipy_wrapper
    elif solver_type == 'petsc4py':
        nonlin_package = petsc4py_wrapper
    
    transient_solve(
        GeometricProperties(
            n_x=n_x,
            n_y=n_y,
            size_x=n_x / 100.,
            size_y=n_y / 100.,
        ),
        PhysicalProperties(
            k=1.e-3,
            rho=1.0,
            c_p=1.0,
        ),
        ConstantMapInitialCondition(
            T=initial_condition,
        ),
        TemperatureBoundaryConditions(
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
        ),
        TimestepProperties(
            delta_t=0.01,
            final_time=0.5,
        ),
        nonlin_package.create(),
        setup_nonlinear_solver=nonlin_package.setup,
    )
