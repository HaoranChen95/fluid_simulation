import numpy as np
from ..file_stream import *


class cell_list:
    def __init__(self, df: "data_file", frame, r_cut_off):
        print("[creating] cell list")
        self.cells = []
        self.cfg = df.cfg[frame]
        self.df_lxyz = df.lxyz
        self.grid_num = df.lxyz // r_cut_off
        self.grid_num += 1 * (self.grid_num == 0)
        self.grid_num = np.int32(self.grid_num)
        gx, gy, gz = self.grid_num
        modified_r_cut_off = df.lxyz / self.grid_num
        for i in range(gx):  # number of x grid
            j_list = []
            for j in range(gy):  # number of y grid
                k_list = []
                for k in range(gz):  # number of z grid
                    k_list.append([])
                j_list.append(k_list)
            self.cells.append(j_list)

        # minium image
        self.cell_xyz = self.cfg - (self.cfg // df.lxyz) * df.lxyz
        self.cell_xyz = np.int32(self.cell_xyz // modified_r_cut_off)
        for i in range(df.Nm):
            gx, gy, gz = self.cell_xyz[i]
            if not np.sum(np.array([gx, gy, gz]) < self.grid_num):
                print(self.cfg)
                raise "d is out of range"
            self.cells[gx][gy][gz].append(i)

    def calc_r_norm(self, i, j):
        r = self.cfg[i, :] - self.cfg[j, :]
        r = ((r + 0.5 * self.df_lxyz) % self.df_lxyz) - 0.5 * self.df_lxyz
        return np.linalg.norm(r)

    def set_ij_list(self):
        self.ij_list = []
        for cell, x, y, z in self.__cell_generator():
            neighbor_c = self.__neighbor_cells(x, y, z)
            for ci, i in enumerate(cell):
                for j in cell[ci + 1 :] + neighbor_c:
                    if i < j:
                        self.ij_list.append((i, j))

    def __cell_generator(self):
        for gx, x_line in enumerate(self.cells):
            for gy, y_line in enumerate(x_line):
                for gz, cell in enumerate(y_line):
                    yield cell, gx, gy, gz

    def __neighbor_cells(self, x, y, z):
        gx, gy, gz = self.grid_num
        neighbor = (
              self.cells[(x + 1) % gx][(y + 1) % gy][(z + 1) % gz]
            + self.cells[(x + 1) % gx][(y + 1) % gy][(z    ) % gz]
            + self.cells[(x + 1) % gx][(y + 1) % gy][(z - 1) % gz]
            + self.cells[(x + 1) % gx][(y    ) % gy][(z + 1) % gz]
            + self.cells[(x + 1) % gx][(y    ) % gy][(z    ) % gz]
            + self.cells[(x + 1) % gx][(y    ) % gy][(z - 1) % gz]
            + self.cells[(x + 1) % gx][(y - 1) % gy][(z + 1) % gz]
            + self.cells[(x + 1) % gx][(y - 1) % gy][(z    ) % gz]
            + self.cells[(x + 1) % gx][(y - 1) % gy][(z - 1) % gz]
            + self.cells[(x    ) % gx][(y + 1) % gy][(z + 1) % gz]
            + self.cells[(x    ) % gx][(y + 1) % gy][(z    ) % gz]
            + self.cells[(x    ) % gx][(y + 1) % gy][(z - 1) % gz]
            + self.cells[(x    ) % gx][(y    ) % gy][(z + 1) % gz]
            + self.cells[(x    ) % gx][(y    ) % gy][(z - 1) % gz]
            + self.cells[(x    ) % gx][(y - 1) % gy][(z + 1) % gz]
            + self.cells[(x    ) % gx][(y - 1) % gy][(z    ) % gz]
            + self.cells[(x    ) % gx][(y - 1) % gy][(z - 1) % gz]
            + self.cells[(x - 1) % gx][(y + 1) % gy][(z + 1) % gz]
            + self.cells[(x - 1) % gx][(y + 1) % gy][(z    ) % gz]
            + self.cells[(x - 1) % gx][(y + 1) % gy][(z - 1) % gz]
            + self.cells[(x - 1) % gx][(y    ) % gy][(z + 1) % gz]
            + self.cells[(x - 1) % gx][(y    ) % gy][(z    ) % gz]
            + self.cells[(x - 1) % gx][(y    ) % gy][(z - 1) % gz]
            + self.cells[(x - 1) % gx][(y - 1) % gy][(z + 1) % gz]
            + self.cells[(x - 1) % gx][(y - 1) % gy][(z    ) % gz]
            + self.cells[(x - 1) % gx][(y - 1) % gy][(z - 1) % gz]
        )
        return neighbor

    # def neighbors_of(self, i):
    #     neighbors = []
    #     cx, cy, cz = self.cell_xyz[i]
    #     for j in self.cells[cx][cy][cz]:
    #         if i != j:
    #             neighbors.append(j)
    #     neighbors += self.__neighbor_cells(cx, cy, cz)
    #     return neighbors
