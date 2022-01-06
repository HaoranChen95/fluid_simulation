import fluid_analyse as fa
import numpy as np

df = fa.data_file()
# fa.write_xyz(df, 1, 100)
fa.total_v_corr(df)
