import numpy as np
from DataAnalyse import _global
# from DataAnalyse.analyse.r_bond import r_bond
from .data_processing import *
from ..algorithm import *
from ..file_stream import *
import matplotlib.pyplot as plt

import multiprocessing as mp


def plot_RDF(bin_hist):
    fig, ax = plt.subplots(dpi=200)

    columns = [r"$r$", r"$p$"]
    label_str = f"""
        Phi = {_global.Phi}
        Pe = {_global.Pe} 
        Nm1p = {_global.Nm1p}
        """.replace("        ", "")
    ax.set_xlabel(columns[0])
    ax.set_ylabel(columns[1])
    ax.plot(
        bin_hist[0][:, 0],
        bin_hist[0][:, 1],
        label=label_str
    )
    # ax.set_yscale("log")
    ax.legend()
    fig.savefig(_global.analyse_dir + f"RDF_{_global.fn_pattern}.jpg")

    plt.close("all")


def RDF():
    hist = mean_counter()
    for frame in _global.frames_selector:
        r_ij = []
        r_norm = []
        _global.cfg_frame = _global.cfg_dset[frame]

        c_list = cell_list(_global.cfg_frame, _global.binrange[0][1])
        

        for cell, x, y, z in c_list.cell_generator():
            neighbor_c = c_list.neighbor_cells(x, y, z)
            for ci, i in enumerate(cell):
                for j in (cell[ci + 1:] + neighbor_c):

                    if i < j and i // _global.Nm1p != j // _global.Nm1p:
                        r_ij.append((i, j))

        # print(r_ij)
        # raise TimeoutError
        with mp.Pool(4) as pool:
            r_norm = pool.starmap(calc_r_ij_norm, r_ij)
        r_norm = np.array(r_norm)
        r_norm = r_norm[r_norm <= _global.binrange[0][1]]
        #! test 
        print("with cell list",len(r_norm))

        # r_ij = []
        # r_norm = []
        # for i in range(_global.Nm):
        #     for j in range(i+1, _global.Nm):
        #         if i // _global.Nm1p != j // _global.Nm1p:
        #             r_ij.append((i,j))
        # with mp.Pool(4) as pool:
        #     r_norm = pool.starmap(calc_r_ij_norm, r_ij)
        # r_norm = np.array(r_norm)
        # r_norm = r_norm[r_norm <= _global.binrange[0][1]]
        # print("without cell list", len(r_norm))
        #! test end
        
        hist += np.histogram(r_norm, bins=_global.get_bins())[0]
        plot_RDF([np.vstack((_global.get_bins_plot(), hist.data)).T])
        frame_time = _global.cfg_data_file["time"][frame]
        write_hist(
            f"RDF_{_global.fn_pattern}.txt",
            [hist.data],
            header=f"frame time: {frame_time}")


def total_RDF():
    print(
        "=" * 80 +
        "\n================ Radius Distribution Function ================\n" +
        "=" * 80)

    if _global.not_load_h5:
        load_h5_cfg_data()

    if not _global.binrange:
        _global.binrange = [(0.5, 6, 0.001)]
        print("set bin range", _global.binrange)
        
    set_cfg_dt()

    RDF()

    _global.binrange = []
