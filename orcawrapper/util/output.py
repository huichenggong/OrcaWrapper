from pathlib import Path
import numpy as np
import pandas as pd

class orca_opt_output:
    def __init__(self, file_name):
        self.file_name = Path(file_name)
        with open(file_name, "r") as f:
            self.lines = f.readlines()

        self.terminated_normally = self.check_terminated_normally()
        self.xyz = self.get_opt_traj()

    def check_terminated_normally(self):
        """
        Check if 'ORCA TERMINATED NORMALLY' is in the last few lines of the output file.
        :return: Ture/False
        """
        for i in [-1, -2, -3]:  # most likely the -2
            if "ORCA TERMINATED NORMALLY" in self.lines[i]:
                return True
        return False

    def get_opt_traj(self):
        """
        Get the optimized trajectory to a list of numpy array (unit in Angstrom)
        read the xyz coordinate after 'CARTESIAN COORDINATES (ANGSTROEM)'
        :return: a list of numpy array
        """
        start = 0
        end = 0
        start_list = []
        end_list = []
        for i, l in enumerate(self.lines):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in l and "GEOMETRY OPTIMIZATION CYCLE" in self.lines[i-3]:
                start = i+2
            elif start and l.strip() == "":
                end = i
                start_list.append(start)
                end_list.append(end)
                start = 0
                end = 0
        xyz_list = []
        for s, e in zip(start_list, end_list):
            xyz_lines = self.lines[s:e]
            xyz = [[float(i) for i in l.split()[1:4]] for l in xyz_lines]
            xyz = np.array(xyz)
            xyz_list.append(xyz)
        return xyz_list


    def get_opt_info(self):
        """
        Get the optimized information
        read the convergence information
        :return:
        """

    def write_traj_to_pdb(self):
        pass

    def write_opt_info_to_csv(self):
        pass

