import numpy as np
import matplotlib.pyplot as plt
from ..file_stream import *


def MSD(df: "data_file"):
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


def plot_MSD(fn_pattern: "str", MSD_data):
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
    time = MSD_data[:, 0]
    ax.plot(time, MSD_data[:, 1])

    ax.legend([label_str], handlelength=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(axis="both", which="both", direction="in")

    ax.grid()
    fig.savefig(ds.dir + f"MSD_{fn.pattern}_grid.jpg")
    plt.close("all")


def total_MSD(df: "data_file"):
    print(
        "=" * 80 + "\n============== Analysing MSD segment ==============\n" + "=" * 80
    )
    df.set_cfg_dt_list()
    MSD_data = np.array(MSD(df)).T
    print(MSD_data.shape)
    ds = data_stream()
    ds.write_data(f"MSD_{df.fn_pattern}.txt", MSD_data.data)

    plot_MSD(df.fn_pattern, MSD_data)
