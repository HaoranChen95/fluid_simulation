from ..file_stream import *
from .color_map import *


def write_xyz(df: "data_file", dt=1, frames=500):
    print("\n ====================== writing cfg xyz ======================")

    df.set_cfg_dt(dt)
    shift_vector = np.array((0.2, 0.2, 0.2))

    cm = color_map()
    with open(f"cfg_{df.fn_pattern}_dt_{dt}.xyz", "a") as f_xzy:
        for frame_i, frame in enumerate(df.frames_index):
            print("[writing] cfg frame", df.time[frame])
            cfg_data = df.cfg[frame, :, :] + shift_vector
            if frame_i == 0:
                color_list = []
                for y in cfg_data[:, 1]:
                    y = round((y % df.lxyz[1]) / df.lxyz[1], 1)
                    color_list.append(cm.blue_grey_red(y))

            f_xzy.write(f"{df.Nm}\n\n")
            cfg_data = cfg_data - (cfg_data // df.lxyz) * df.lxyz

            for i in range(df.Nm):
                write_str = ""
                for val in color_list[i]:
                    write_str += f"{val:.4f}, "
                for val in cfg_data[i]:
                    write_str += f"{val:.2f}, "
                write_str += "\n"
                f_xzy.write(write_str)

            if frame_i >= frames:
                break
