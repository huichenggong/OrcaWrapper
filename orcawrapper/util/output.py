from pathlib import Path
import numpy as np
import pandas as pd
import re

def total_run_time(line):
    """
    Given the orca final output line, return the total run time in seconds.
    :param line:
    :return: float, time in seconds
    """
    words = line.split()

    days = float(words[3])
    hours = float(words[5])
    minutes = float(words[7])
    seconds = float(words[9])
    msec = float(words[11])
    time = (((days * 24) + hours) * 60 + minutes) * 60 + seconds + msec / 1000
    return time

class orca_output:
    def __init__(self, file_name):
        self.file_name = Path(file_name)
        with open(file_name, "r") as f:
            self.lines = f.readlines()
        self.content = "".join(self.lines)
        self.time = None
        self.terminated_normally = self.check_terminated_normally()

        self.basic_info = {}
        self.get_basic_info()

    def check_terminated_normally(self):
        """
        Check if 'ORCA TERMINATED NORMALLY' is in the last few lines of the output file.
        :return: Ture/False
        """
        for i in [-1, -2, -3, -4]:  # most likely the -2
            if "ORCA TERMINATED NORMALLY" in self.lines[i]:
                self.time = total_run_time(self.lines[i+1])
                return True
        return False

    def get_basic_info(self):
        """
        Get the basic information of the output file
        :return: a dictionary
        """
        self.basic_info = {}
        match = re.search(r'Number of atoms\s*\.*\s*(\d+)', self.content)
        if match:
            self.basic_info["n_atoms"] = int(match.group(1))
        match = re.search(r'Number of basis functions\s*\.*\s*(\d+)', self.content)
        if match:
            self.basic_info["n_basis_functions"] = int(match.group(1))



class orca_opt_output(orca_output):
    def __init__(self, file_name):
        orca_output.__init__(self, file_name)

        self.xyz = self.get_opt_traj()
        self.convergence_info_df = self.get_opt_info()

    def check_terminated_normally(self):
        """
        Check if 'ORCA TERMINATED NORMALLY' is in the last few lines of the output file.
        :return: Ture/False
        """
        for i in [-1, -2, -3, -4]:  # most likely the -2
            if "ORCA TERMINATED NORMALLY" in self.lines[i]:
                self.time = total_run_time(self.lines[i+1])
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
        :return: a pandas DataFrame
        """
        Energy_delta = []
        RMS_gradient = []
        MAX_gradient = []
        RMS_step = []
        MAX_step = []
        for i, l in enumerate(self.lines):
            if "Geometry convergence" in l:
                # check the incoming i+3 - i+7 lines
                e = 0
                rms_grad = 0
                max_grad = 0
                rms_step = 0
                max_step = 0
                for line in self.lines[i+3:i+8]:
                    if "Energy change" in line:
                        e = float(line.split()[2])
                    elif "RMS gradient" in line:
                        rms_grad = float(line.split()[2])
                    elif "MAX gradient" in line:
                        max_grad = float(line.split()[2])
                    elif "RMS step" in line:
                        rms_step = float(line.split()[2])
                    elif "MAX step" in line:
                        max_step = float(line.split()[2])
                Energy_delta.append(e)
                RMS_gradient.append(rms_grad)
                MAX_gradient.append(max_grad)
                RMS_step.append(rms_step)
                MAX_step.append(max_step)
        df = pd.DataFrame({"Energy_delta": Energy_delta,
                           "RMS_gradient": RMS_gradient,
                           "MAX_gradient": MAX_gradient,
                           "RMS_step": RMS_step,
                           "MAX_step": MAX_step})
        return df


    def write_traj_to_pdb(self):
        pass

    def write_opt_info_to_csv(self):
        pass

