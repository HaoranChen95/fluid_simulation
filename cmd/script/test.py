import fluid_analyse as fa
import numpy as np

df = fa.data_file(lxyz=(40, 40, 40))
# fa.write_xyz(df, 1, 100)
fa.total_MSD(df)

print(df.time[-10:])
print(len(df.time))
df.set_cfg_dt(1)
print(df.time[df.frames_index][-10:])