import os
from .file_name import *


class file_list:
    def __init__(self, head, suffix, dir="./", sort_by="aV") -> None:
        self.list = [
            dir + x
            for x in os.listdir(dir)
            if (x.find(head) == 0 and x.find(suffix) > 0)
        ]
        self.list.sort(key=lambda x: find_param_float(sort_by, x))

    def print(self):
        for fn in self.list:
            print(fn)
