import numpy as np
from time import time
from solver.properties import GeometricProperties, PhysicalProperties,\
    ConstantInitialCondition, TimestepProperties, PrescribedTemperatureBoundaryCondition,\
    TemperatureBoundaryConditions, PrescribedFlowBoundaryCondition, FromFileInitialCondition
from solver.nonlinear_solver import solve


def plot_animated_results(result_list):
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    N_FRAMES = len(result_list)
    x, y = result_list[0].shape
    x = np.arange(0, x)
    y = np.arange(0, y)

    f = plt.figure()
    img = plt.imshow(result_list[0], vmin=np.min(result_list), vmax=np.max(result_list))

    def animate(i):
        img.set_data(result_list[i])

    anim = animation.FuncAnimation(
        f,
        animate,
        frames=N_FRAMES,
        interval=60,
        repeat_delay=1000
    )

    plt.show()


class InMemoryResultsHandler():

    def __init__(self):
        self._s = time()
        self._i = 0
        self.result_list = []

    def __call__(self, time_, result):
        print("Saving time[%s] %2.5f..." % (self._i, time_,), end='')
        self.result_list.append(np.array(result))
        print(" Done. (%s)" % (time() - self._s,))
        self._s = time()
        self._i += 1


def test_solver():
    result_handler = InMemoryResultsHandler()
    n_x = 200
    n_y = 50

    solve(
        GeometricProperties(
            n_x=n_x,
            n_y=n_y,
            size_x=1.0,
            size_y=0.1,
        ),
        PhysicalProperties(
            k=10.0,
            rho=1.0,
            c_p=1.0,
        ),
        ConstantInitialCondition(
            T=0.0,
            n_x=n_x,
            n_y=n_y,
        ),
        TemperatureBoundaryConditions(
            PrescribedTemperatureBoundaryCondition(lambda t: 0.0),
#             PrescribedFlowBoundaryCondition(lambda t: 0.0),
            PrescribedTemperatureBoundaryCondition(lambda t: np.sin(t*10)),
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
            PrescribedFlowBoundaryCondition(lambda t: 0.0),
        ),
        TimestepProperties(
            delta_t=0.01,
            final_time=2.0,
        ),
        on_timestep_callback=result_handler,
    )

    plot_animated_results(result_handler.result_list)
