import fluid_analyse as fa
import numpy as np

df = fa.data_file()
fa.write_xyz(df, 0.01, 1000)
# fa.total_v_dist(df)
