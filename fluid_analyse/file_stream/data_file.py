import numpy as np
import h5py
import os

from .file_name import *
from .file_list import *


class data_file:
    def __init__(self, dir="./", lxyz=(20, 20, 20)) -> None:
        fl = file_list("cfg_data_", ".txt", dir)
        if fl.list:
            self.convert_txt_to_hdf5()
        fl = file_list("cfg_data_", ".h5", dir)
        fn = file_name(fl.list[0])
        self.Nm = fn.Nm
        self.phi = fn.phi
        self.kT = fn.kT
        self.gamma = fn.gamma
        self.fn_pattern = fn.pattern
        self.lxyz = np.array(lxyz)
        self.file = h5py.File(fn.str, "r")
        self.cfg = self.file["cfg"]
        self.vel = self.file["vel"]
        self.time = np.array(self.file["time"])

    def convert_txt_to_hdf5(self):
        print("calling convert_txt_to_hdf5")
        fl = file_list("cfg_data", ".txt")
        if not fl.list:
            print("no data to convert")
            return 0
        fn = file_name(fl.list[0])

        with h5py.File(f"cfg_data_{fn.pattern}.h5", "w") as f:

            print("[loading]", fl.list[0])
            data = np.loadtxt(fl.list[0])
            print(data.shape)
            print(data[:5, :3])
            n_xyz = len(data[0])
            data = np.reshape(data, (-1, fn.Nm, n_xyz))
            print(data.shape)
            print(data[0, :5, :3])
            dset_names = ["cfg", "vel"]
            for i in range(n_xyz // 3):
                f.create_dataset(
                    dset_names[i],
                    data=data[:, :, i * 3 : (i + 1) * 3],
                    compression="gzip",
                    chunks=(500, 1, 3),
                    maxshape=(None, None, 3),
                )

            for fn_str in fl.list[1:]:
                fn = file_name(fn_str)
                print("[loading]", fn_str)
                data = np.loadtxt(fn_str)
                print(data.shape)
                print(data[:-5, :3])
                data = np.reshape(data, (-1, fn.Nm, n_xyz))
                print(data.shape)
                print(data[-1, :-5, :3])
                for i in range(n_xyz // 3):
                    h5_dset_append(f[dset_names[i]], data[:, :, i * 3 : (i + 1) * 3])
                print(f["cfg"].shape)
                print(f["cfg"][-1, :-5, :3])

            time_data = self.__load_time()
            assert len(time_data) == len(f["cfg"])
            f.create_dataset("time", data=time_data, compression="gzip")

            print("[removing] cfg_data_*.txt")
            for fn_str in file_list("cfg_data", ".txt", sort_by="aV").list:
                os.remove(fn_str)
    
    def set_cfg_dt(self, dt=1):
        print("[setting] cfg_dt with frames at time: ")
        time = 0
        self.frames_index = []
        for i, t in enumerate(self.time):
            if (t - time) ** 2 < 1e-8:
                self.frames_index.append(i)
                time += dt
        print(self.time[self.frames_index][:10], "...")
        print(self.time[self.frames_index][-10:])

    def set_cfg_dt_list(self, periodic=10, amount=100):
        print("[getting] dt series with periodic", periodic, "and amount", amount)
        print(self.time[:10])
        print(self.time[-10:])

        self.frames_index_list = []

        start_indexes = np.searchsorted(
            self.time, list(range(0, periodic * amount, periodic))
        )
        print(self.time[start_indexes])

        t_intersection = self.time
        half_dt = 0.5 * t_intersection[1] - t_intersection[0]
        for start in start_indexes[1:]:
            print(len(t_intersection))
            t_intersection_index = np.searchsorted(
                t_intersection,
                self.time[start:]
                - self.time[start]
                - half_dt,
            )
            t_intersection = t_intersection[t_intersection_index]
            increase = t_intersection[1:] > t_intersection[:-1]
            while not np.any(increase):
                t_intersection = np.append(t_intersection[:-1][increase], t_intersection[-1])
                increase = t_intersection[1:] > t_intersection[:-1]
        print(len(t_intersection))

        for start_index in start_indexes:
            time = self.time[start_index] + t_intersection - half_dt
            index = start_index
            time_indexes = []
            for t in time:
                while t > self.time[index]:
                    index += 1
                time_indexes.append(index)
            self.frames_index_list.append(time_indexes)

        # print the time value of index
        for time_indexes in self.frames_index_list:
            print(self.time[time_indexes][:10])
            print(self.time[time_indexes][-10:])

    def __load_time(self):
        print("[loading] the time file")
        fl = file_list("cfg_time_", ".txt", sort_by="aV")
        if fl:
            data = []
            for fn in fl.list:
                data = np.hstack((data, np.loadtxt(fn)))
            return data


def h5_dset_append(dset, other_dset):
    dset.resize(dset.shape[0] + other_dset.shape[0], axis=0)
    dset[-other_dset.shape[0] :] = other_dset
