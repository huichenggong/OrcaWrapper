import os
from pathlib import Path
import numpy as np

class orca_input:
    def __init__(self, file_name):
        self.file_name = Path(file_name)
        with open(file_name, "r") as f:
            self.lines = f.readlines()
        b1, ind_tuple = self.check_xyz()
        if b1:
            self.xyz = self.get_xyz(ind_tuple)
        else:
            self.xyz = None


    def __str__(self):
        return f"ORCA input file : {self.file_name.__str__()}"

    def __repr__(self):
        return self.file_name.__str__()

    def check_xyz(self):
        """
        Check if xyz (coordinate) is provided in the input file.
        coordinate can be provided in xyz section or separate file
        :return:
        """
        xyz_start = 0
        xyz_end = 0
        is_xyz = False
        for i, l in enumerate(self.lines):
            if ["*", "xyz"] == l.split()[:2]:
                is_xyz = True
                xyz_start = i
            elif is_xyz and l.split()[0] == "*":
                xyz_end = i
                return True, (xyz_start, xyz_end)
        return False, 0

    def get_xyz(self, ind_tuple):
        """
        read xyz from input file
        :param ind_tuple: tuple of start and end index of xyz section
        :return: xyz coordinate in numpy array
        """
        xyz_lines = self.lines[ind_tuple[0]+1: ind_tuple[1]]
        xyz = [[float(i) for i in l.split()[1:4]] for l in xyz_lines]
        xyz = np.array(xyz)
        return xyz
