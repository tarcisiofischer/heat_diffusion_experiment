from time import time
import numpy as np

class AnimationPlotterResultsHandler(object):

    def __init__(self, verbose=False):
        self._s = time()
        self._i = 0
        self.result_list = []
        self.verbose = verbose


    def __call__(self, time_, result):
        if self.verbose:
            print("Saving time[%s] %2.5f..." % (self._i, time_,), end='')

        self.result_list.append(np.array(result))

        if self.verbose:
            print(" Done. (%s)" % (time() - self._s,))

        self._s = time()
        self._i += 1


    def plot(self):
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
    
        N_FRAMES = len(self.result_list)
        x, y = self.result_list[0].shape
        x = np.arange(0, x)
        y = np.arange(0, y)
    
        f = plt.figure()
        vmin = np.min(self.result_list)
        vmax = np.max(self.result_list)
        img = plt.imshow(self.result_list[0], vmin=vmin, vmax=vmax)
    
        def animate(i):
            img.set_data(self.result_list[i])
    
        # Although the _anim variable is unused, it must be kept alive in
        # order to the mechanism to work.
        _anim = animation.FuncAnimation(
            f,
            animate,
            frames=N_FRAMES,
            interval=60,
            repeat_delay=1000
        )
    
        plt.show()
