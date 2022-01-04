import numpy as np
from .. import _global


def minimum_image(r):
    print(_global.lx, _global.ly, _global.lz)


def calc_r_ij_norm(i, j):
    r = (
        _global.cfg_frame[i, :]
        -
        _global.cfg_frame[j, :]
    )
    r = (
        ((r + 0.5 * _global.lxyz)
            % _global.lxyz)
        - 0.5 * _global.lxyz
    )
    return np.linalg.norm(r)


class cell_list:
    def __init__(self, cfg_data, r_cut_off):
        print("[creating] cell list")
        self.cells = []
        self.grid_num = (_global.lxyz // r_cut_off)
        self.grid_num += 1 * (self.grid_num == 0)
        self.grid_num = np.int32(self.grid_num)
        gx, gy, gz = self.grid_num
        modified_r_cut_off = _global.lxyz / self.grid_num
        for i in range(gx):  # number of x grid
            j_list = []
            for j in range(gy):  # number of y grid
                k_list = []
                for k in range(gz):  # number of z grid
                    k_list.append([])
                j_list.append(k_list)
            self.cells.append(j_list)

        # minium image
        self.cell_xyz = cfg_data - (
            cfg_data // _global.lxyz) * _global.lxyz
        self.cell_xyz = np.int32(self.cell_xyz // modified_r_cut_off)
        for i in range(_global.Nm):
            gx, gy, gz = self.cell_xyz[i]
            if not np.sum(np.array([gx, gy, gz]) < self.grid_num):
                print(cfg_data)
                raise "d is out of range"
            self.cells[gx][gy][gz].append(i)

    def cell_generator(self):

        for gx, x_line in enumerate(self.cells):
            for gy, y_line in enumerate(x_line):
                for gz, cell in enumerate(y_line):
                    yield cell, gx, gy, gz

    def half_neighbor_cells(self, x, y, z):
        gx, gy, gz = self.grid_num
        if gz == 1:
            neighbor = (
                self.cells[(x + 1) % gx][(y + 1) % gy][0] +
                self.cells[x][(y + 1) % gy][0] +
                self.cells[(x - 1) % gx][(y + 1) % gy][0] +
                self.cells[(x + 1) % gx][y][0]
            )
        elif gz == 3:
            neighbor = (
                self.cells[(x + 1) % gx][(y + 1) % gy][z] +
                self.cells[x][(y + 1) % gy][z] +
                self.cells[(x - 1) % gx][(y + 1) % gy][z] +
                self.cells[(x + 1) % gx][y][z] +
                self.cells[(x + 1) % gx][y][z - 1] +
                self.cells[(x) % gx][y][z - 1] +
                self.cells[(x - 1) % gx][y][z - 1] +
                self.cells[(x + 1) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x - 1) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x + 1) % gx][(y - 1) % gy][z - 1] +
                self.cells[(x) % gx][(y - 1) % gy][z - 1] +
                self.cells[(x - 1) % gx][(y - 1) % gy][z - 1]
            )
        return neighbor

    def neighbor_cells(self, x, y, z):
        gx, gy, gz = self.grid_num
        if gz == 1:
            neighbor = (
                self.cells[(x - 1) % gx][(y + 1) % gy][0] +
                self.cells[x][(y + 1) % gy][0] +
                self.cells[(x + 1) % gx][(y + 1) % gy][0] +
                self.cells[(x - 1) % gx][(y - 1) % gy][0] +
                self.cells[x][(y - 1) % gy][0] +
                self.cells[(x + 1) % gx][(y - 1) % gy][0] +
                self.cells[(x + 1) % gx][y][0] +
                self.cells[(x - 1) % gx][y][0]
            )
        elif gz == 3:
            print("just half finished")
            raise TimeoutError
            neighbor = (
                self.cells[(x + 1) % gx][(y + 1) % gy][z] +
                self.cells[x][(y + 1) % gy][z] +
                self.cells[(x - 1) % gx][(y + 1) % gy][z] +
                self.cells[(x + 1) % gx][y][z] +
                self.cells[(x + 1) % gx][y][z - 1] +
                self.cells[(x) % gx][y][z - 1] +
                self.cells[(x - 1) % gx][y][z - 1] +
                self.cells[(x + 1) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x - 1) % gx][(y + 1) % gy][z - 1] +
                self.cells[(x + 1) % gx][(y - 1) % gy][z - 1] +
                self.cells[(x) % gx][(y - 1) % gy][z - 1] +
                self.cells[(x - 1) % gx][(y - 1) % gy][z - 1]
            )
        return neighbor

    def neighbors_of(self, i):
        neighbors = []
        cx, cy, cz = self.cell_xyz[i]
        for j in self.cells[cx][cy][cz]:
            if i != j:
                neighbors.append(j)
        neighbors += self.neighbor_cells(cx,cy,cz)
        return neighbors
