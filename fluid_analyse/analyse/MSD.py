from DataAnalyse import _global
from ..file_stream import *
from ..plot import *
from ..algorithm import *
import numpy as np
from .data_processing import *

def MSD_segment(index_list):
    time_list = _global.cfg_data_file["time"][index_list[1:]] - _global.cfg_data_file["time"][index_list[0]]
    r_sq = np.sum((_global.cfg_dset[index_list[1:], :, :] -
                   _global.cfg_dset[index_list[0], :, :][np.newaxis, :, :]
                   )**2,
                  axis=2)
    MSD = np.mean(r_sq, axis=1)
    if _global.Nm1p == 1:
        return time_list, MSD
    MSD_cm = np.mean(
        np.sum((
            _global.cfg_cm_dset[index_list[1:], :, :] -
            _global.cfg_cm_dset[index_list[0], :, :][np.newaxis, :, :])**2,
            axis=2),
        axis=1)

    return time_list, MSD, MSD_cm


def total_MSD_segment():
    # _global.MSD_dt = 1
    # _global.MSD_t_end = 500
    print("=" * 80 +
          "\n============== Analysing MSD segment ==============\n" + "=" * 80)
    if _global.not_load_h5:
        load_h5_cfg_data()

    if _global.Nm1p != 1:
        get_center_of_mass()
    # print(np.mean(_global.cfg_cm_dset[:,:,:50],axis=0))

    print("[calculate] MSD segment")
    indexes_list = get_cfg_dt_list()
    MSD_segment_data = mean_counter()
    for index_list in indexes_list:
        MSD_segment_data += np.array(MSD_segment(index_list)).T
    print(MSD_segment_data.data.shape)
    write_data(f"MSD_segment_{_global.fn_pattern}.txt", MSD_segment_data.data)

    r_g_square = 0
    if _global.Nm1p != 1:
        r_g_square = read_header(
            file_list_generator("r_gyration", ".txt", "./analyse/")[0]
        )
        if r_g_square:
            r_g_square = r_g_square[1][0]

    plot_MSD_segment(MSD_segment_data.data, r_g_square)
