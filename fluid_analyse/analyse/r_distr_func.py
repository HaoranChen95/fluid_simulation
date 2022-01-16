import numpy as np
from ..file_stream import *
from .cell_list import *
import matplotlib.pyplot as plt

import multiprocessing as mp


def plot_RDF(fn_pattern: "str", data):
    ds = data_stream()
    fn = file_name(fn_pattern)
    fig, ax = plt.subplots(dpi=200)

    columns = [r"$r$", r"$p$"]
    label_str = f"""
        $\phi$ = {fn.phi}
        $\gamma$ = {fn.gamma}
        """.replace(
        "        ", ""
    )
    ax.set_xlabel(columns[0])
    ax.set_ylabel(columns[1])
    ax.plot(data[:, 0], data[:, 1], label=label_str)
    ax.legend()
    fig.savefig(ds.dir + f"RDF_{fn.pattern}.jpg")

    plt.close("all")


def RDF(df: "data_file"):
    hist = mean_counter()
    br = binrange(0.5, 6, 0.001)
    for frame in df.frames_index:
        r_norm = []

        c_list = cell_list(df, frame, br.end)
        c_list.set_ij_list()

        with mp.Pool(4) as pool:
            r_norm = pool.starmap(c_list.calc_r_norm, c_list.ij_list)
        r_norm = np.array(r_norm)
        r_norm = r_norm[r_norm <= br.end]

        # #! test start
        # print("with cell list", len(r_norm))
        # r_ij = []
        # r_norm = []
        # for i in range(df.Nm):
        #     for j in range(i + 1, df.Nm):
        #         r_ij.append((i, j))
        # with mp.Pool(4) as pool:
        #     r_norm = pool.starmap(c_list.calc_r_norm, r_ij)
        # r_norm = np.array(r_norm)
        # r_norm = r_norm[r_norm <= br.end]
        # print("without cell list", len(r_norm))
        # #! test end

        hist += np.histogram(r_norm, bins=br.get_bins())[0]
        data = np.vstack((br.get_bins_plot(), hist.data)).T
        plot_RDF(df.fn_pattern, data)
        frame_time = df.time[frame]
        ds = data_stream()
        ds.write_data(
            f"RDF_{df.fn_pattern}.txt",
            data,
            header=f"frame time: {frame_time}",
        )


def total_RDF(df: "data_file"):
    print(
        "=" * 80
        + "\n================ Radius Distribution Function ================\n"
        + "=" * 80
    )

    df.set_cfg_dt(1)

    RDF(df)
