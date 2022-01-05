import numpy as np
import matplotlib.pyplot as plt
from ..file_stream import *


def v_dist(df: "data_file"):
    resalt = mean_counter()
    time_list = (
        df.time[df.frames_index_list[0][1:]] - df.time[df.frames_index_list[0][0]]
    )
    for i_list in df.frames_index_list:
        print(df.time[i_list[0]])
        data = np.mean(
            np.sum((df.cfg[i_list[1:], :, :] - df.cfg[i_list[0], :, :]) ** 2, axis=2),
            axis=1,
        )
        resalt += np.array(data)

    return time_list, resalt.data


def plot_v_dist(fn_pattern: "str", v_dist_data):
    ds = data_stream()
    fn = file_name(fn_pattern)

    columns = [r"$t$", r"$\langle \vec{r}(t)^2 \rangle$"]
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
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(axis="both", which="both", direction="in")
    fig.savefig(ds.dir + f"v_dist_{fn.pattern}.jpg")

    ax.grid()
    fig.savefig(ds.dir + f"v_dist_{fn.pattern}_grid.jpg")
    plt.close("all")


def total_v_dist(df: "data_file"):
    print(
        "=" * 80 + "\n============== Analysing v_dist segment ==============\n" + "=" * 80
    )
    df.set_cfg_dt_list()
    v_dist_data = np.array(v_dist(df)).T
    print(v_dist_data.shape)
    ds = data_stream()
    ds.write_data(f"v_dist_{df.fn_pattern}.txt", v_dist_data.data)

    plot_v_dist(df.fn_pattern, v_dist_data)
