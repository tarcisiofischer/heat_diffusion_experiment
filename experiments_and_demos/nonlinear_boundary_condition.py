from solver.output_handlers.animation_plotter import AnimationPlotterResultsHandler
import numpy as np
from solver.nonlinear_solver_wrappers import scipy_wrapper,\
    petsc4py_wrapper
from solver.time_stepper import transient_solve
from solver.properties import GeometricProperties, PhysicalProperties,\
    ConstantMapInitialCondition, TemperatureBoundaryConditions,\
    PrescribedFlowBoundaryCondition, TimestepProperties,\
    PrescribedTemperatureBoundaryCondition

solver_type = 'scipy'

result_handler = AnimationPlotterResultsHandler(verbose=True)
n_x = 100
n_y = 100

initial_condition = np.zeros(shape=(n_y, n_x))

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
        k=10.0,
        rho=10.0,
        c_p=1.0,
    ),
    ConstantMapInitialCondition(
        T=initial_condition,
    ),
    TemperatureBoundaryConditions(
        east_bc=PrescribedTemperatureBoundaryCondition(lambda t: np.sin(10*t)),
        west_bc=PrescribedTemperatureBoundaryCondition(lambda t: np.cos(10*t)),
        north_bc=PrescribedFlowBoundaryCondition(lambda t: 0.0),
        south_bc=PrescribedFlowBoundaryCondition(lambda t: 0.0),
    ),
    TimestepProperties(
        delta_t=0.01,
        final_time=1.0,
    ),
    nonlin_package.create(),
    setup_nonlinear_solver=nonlin_package.setup,
    on_timestep_callback=result_handler,
)

result_handler.plot()
