import numpy as np
import matplotlib.pyplot as plt
from ..file_stream import *


def v_dist(df: "data_file"):
    br = binrange(-10, 10, 0.01)
    data = np.histogram(df.vel[df.frames_index], br.get_bins(),density=True)[0]
    print(data.shape, data[::100])

    return br.get_bins_plot(), data


def plot_v_dist(fn_pattern: "str", v_dist_data):
    ds = data_stream()
    fn = file_name(fn_pattern)

    columns = [r"$v$", r"$\rho$"]
    label_str = f"""
                Density = {fn.phi} 
                $\gamma$ = {fn.gamma} 
                """.replace(
        "        ", ""
    )

    fig, ax = plt.subplots(dpi=200)

    ax.set_xlabel(columns[0])
    ax.set_ylabel(columns[1])
    time = v_dist_data[:, 0]
    ax.plot(time, v_dist_data[:, 1])

    ax.legend([label_str], handlelength=0)
    ax.tick_params(axis="both", which="both", direction="in")
    ax.grid()
    fig.savefig(ds.dir + f"v_dist_{fn.pattern}.jpg")
    plt.close("all")


def total_v_dist(df: "data_file"):
    print(
        "=" * 80 + "\n============== Analysing v_dist segment ==============\n" + "=" * 80
    )
    df.set_cfg_dt(1)
    v_dist_data = np.array(v_dist(df)).T
    print(v_dist_data.shape)
    ds = data_stream()
    ds.write_data(f"v_dist_{df.fn_pattern}.txt", v_dist_data.data)

    plot_v_dist(df.fn_pattern, v_dist_data)
