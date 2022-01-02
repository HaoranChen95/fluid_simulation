import numpy as np
import h5py

from .file_name import *

class data_file:
    def __init__(self, fn:"file_name") -> None:
        self.fn = fn
        
