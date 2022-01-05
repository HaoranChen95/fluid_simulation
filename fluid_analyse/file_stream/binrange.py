import numpy as np


class binrange:
    def __init__(self, start, end, step) -> None:
        self.start = start
        self.end = end
        self.step = step

    def get_bins(self):
        return np.arange(
            self.start,
            self.end + self.step,
            self.step,
        )

    def get_bins_plot(self):
        return (
            np.arange(
                self.start,
                self.end + self.step,
                self.step,
            )[:-1]
            + 0.5 * self.step
        )
