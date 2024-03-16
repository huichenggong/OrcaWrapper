import os
from pathlib import Path
import numpy as np

class orca_input:
    def __init__(self, file_name):
        self.file_name = Path(file_name)
        with open(file_name, "r") as f:
            self.lines = f.readlines()
        self.has_xyz, self.xyz_index_tuple = self.check_xyz()
        if self.has_xyz:
            self.xyz = self.get_xyz(self.xyz_index_tuple)
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

    def update_xyz(self, xyz):
        """
        update the xyz coordinate in the input file
        :param xyz: 3xN np.array, new xyz coordinate
        :return: None
        """
        if self.xyz.shape != xyz.shape:
            raise ValueError("New xyz shape does not match with the old one.")
        else:
            self.xyz = xyz

    def write(self, file_name, xyz_format="15.12f"):
        """
        write the orca input file. We will use the same format as the original file. and write new xyz.
        :param file_name: str, file name
        :return: None
        """
        new_lines = []
        # lines before xyz can be copied as it is
        for l in self.lines[:self.xyz_index_tuple[0]+1]:
            new_lines.append(l)

        # write xyz
        for l_old, xyz in zip(self.lines[self.xyz_index_tuple[0]+1: self.xyz_index_tuple[1]], self.xyz):
            words = l_old.split()
            if len(words) == 4:
                l_new = l_old.split()[0] + " " + " ".join([f"{i:{xyz_format}}" for i in xyz]) + "\n"
            elif len(words) > 4:
                l_new = l_old.split()[0] + " " + " ".join([f"{i:{xyz_format}}" for i in xyz]) + " " + " ".join(words[4:]) + "\n"
            else:
                raise ValueError("xyz line has less than 4 words.")
            new_lines.append(l_new)

        # lines after xyz can be copied as it is
        for l in self.lines[self.xyz_index_tuple[1]:]:
            new_lines.append(l)

        with open(file_name, "w") as f:
            f.writelines(new_lines)


