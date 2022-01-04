import os
import numpy as np

class data_stream:
    def __init__(self, dir="analyse/") -> None:
        self.dir = dir
        if not os.path.isdir(dir):
            os.mkdir(dir)
            print(f"create directory {dir}")

    def write_data(self, fn, data, header=""):
        print("\n ====================== writing normal data ======================")
        data = np.array(data)
        print(
            f"write data with shape {data.shape} to \n {self.dir + fn}"
        )
        print(f"whith the header:\n{header}")
        np.savetxt(self.dir + fn, data, header=header)
