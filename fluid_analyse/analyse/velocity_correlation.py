import numpy as np
import matplotlib.pyplot as plt
from ..file_stream import *


def v_corr(df: "data_file"):
    resalt = mean_counter()
    time_list = df.time[df.frames_index_list[0]] - df.time[df.frames_index_list[0][0]]
    for i_list in df.frames_index_list:
        print(df.time[i_list[0]])

        resalt += np.mean(
            np.sum((df.vel[i_list, :, :] * df.vel[i_list[0], :, :]), axis=2),
            axis=1,
        )

        data = np.vstack((time_list, resalt.data)).T
        plot_v_corr(df.fn_pattern, data)
        frame_time = df.time[i_list[0]]
        ds = data_stream()
        ds.write_data(
            f"v_corr_{df.fn_pattern}.txt",
            data,
            header=f"frame time: {frame_time}",
        )


def plot_v_corr(fn_pattern: "str", v_corr_data):
    ds = data_stream()
    fn = file_name(fn_pattern)

    columns = [
        r"$t$",
        r"$\langle \vec{v}(0)\cdot\vec{v}(t) \rangle / \langle \vec{v}(0)^2 \rangle$",
    ]
    label_str = f"""
                Density = {fn.phi} 
                $\gamma$ = {fn.gamma} 
                """.replace(
        "        ", ""
    )

    fig, ax = plt.subplots(dpi=200)

    ax.set_xlabel(columns[0])
    ax.set_ylabel(columns[1])
    time = v_corr_data[:, 0]
    ax.plot(time, v_corr_data[:, 1] / v_corr_data[0, 1])

    ax.legend([label_str], handlelength=0)
    ax.tick_params(axis="both", which="both", direction="in")
    
    ax.grid()
    fig.savefig(ds.dir + f"v_corr_{fn.pattern}_grid.jpg")
    plt.close("all")


def total_v_corr(df: "data_file"):
    print(
        "=" * 80
        + "\n============== Analysing v_corr segment ==============\n"
        + "=" * 80
    )
    df.set_cfg_dt_list()
    v_corr(df)
