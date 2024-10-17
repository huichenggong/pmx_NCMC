#!/usr/bin/env python3
import argparse
import sys
import os
import shutil
from pathlib import Path
import time
from functools import wraps
import logging
import tempfile
import subprocess

import numpy as np
import pandas as pd

import pmxNCMC
from pmxNCMC import util


fun_exe_times = {
    'run_eq_grompp': [],
    'run_eq_mdrun' : [],
    'run_ti_grompp': [],
    'run_ti_mdrun' : [],
}
def time_function(func_name):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            duration = end_time - start_time
            fun_exe_times[func_name].append(duration)
            return result
        return wrapper
    return decorator


def logging_ave_time_cycle(ave_time_cycle, time_cycle):
    """
    :param ave_time_cycle: time in seconds
        average time per cycle
    :param time_cycle: time in seconds
        time for this cycle
    :return:
    """
    if ave_time_cycle > 3600:
        ave_time_cycle_str = f"{int(ave_time_cycle // 3600)} h {int((ave_time_cycle % 3600) // 60)} min {(ave_time_cycle % 60):.0f} s"
    elif ave_time_cycle > 60:
        ave_time_cycle_str = f"{int(ave_time_cycle // 60)} min {(ave_time_cycle % 60):.1f} s"
    else:
        ave_time_cycle_str = f"{ave_time_cycle:.2f} s"
    logging.info(f"Time of this cycle : {time_cycle:.2f} s. Average time per cycle: {ave_time_cycle_str}")

class PMX_MDRUN_RE:
    def __init__(self, top, csv, mdp_folder, folder_start, MDRUN, GROMPP, tmp_folder, env, min_output,
                 debug=False):
        """
        :param top: pathlib.Path, gromacs topology file
        :param csv: pathlib.Path, csv file to save work values
        :param mdp_folder: pathlib.Path, folder that contains eq0.mdp, eq1.mdp, ti0.mdp, ti1.mdp
        :param folder_start: pathlib.Path, folder to start the simulation.
        :param MDRUN: str, command for mdrun, we will use multidir, MPI is required.
            For example "mpirun -np 2 gmx_mpi mdrun"
        :param GROMPP: str, command for grompp, with additional flags.
            For example "gmx_threads_AVX_256 grompp -maxwarn 1"
        :param tmp_folder: pathlib.Path, temporary folder where the TI job will run.
        :param env: dict, environment variables
        :param min_output: bool, set to True will redirect all gromacs stdout/stderr to /dev/null
        :param debug: bool, set to True will save ti files

        restart a simulation from 000000 is not possible. It will be identified as a new start.
        """
        self.top           = top
        self.csv           = csv
        self.mdp_folder    = mdp_folder
        self.folder_start  = folder_start
        self.MDRUN         = MDRUN
        self.GROMPP        = GROMPP
        self.tmp_folder    = tmp_folder
        self.env           = env
        self.min_output    = min_output
        self.debug         = debug

        self.safe_flag = False
        self.s0 = [None, None] # starting cpt/gro files for the next eq in lambda 0
        self.s1 = [None, None] # starting cpt/gro files for the next eq in lambda 1

        self.current_cycle = None
        self.ref_t = None
        self.kBT = None
        self.current_folder = None

    def set_current_folder(self):
        """
        set self.current_folder based on self.current_cycle
        At the moment 06d
        This will not create the folder.
        :return: None
        """
        self.current_folder = Path(f"{self.current_cycle:06d}")

    def prepare_current_folder(self):
        """
        prepare the current_folder. Make sure 0 and 1 folders exist.
        :return:
        """
        (self.current_folder / "0").mkdir(parents=True, exist_ok=True)
        (self.current_folder / "1").mkdir(parents=True, exist_ok=True)

    def prepare_scratch_folder(self):
        """
        prepare the scratch folder. Make sure 0 and 1 folders exist.
        :return:
        """
        (self.tmp_folder / "0").mkdir(parents=True, exist_ok=True)
        (self.tmp_folder / "1").mkdir(parents=True, exist_ok=True)
        # make sure both folders are empty
        for i in range(2):
            tmp_ = self.tmp_folder / str(i)
            for f in tmp_.iterdir():
                if f.is_file():
                    f.unlink()
                else:
                    shutil.rmtree(f)

    def safety_check(self):
        """
        Temperature should be the same in all 4 mdp files.
        nstep * delta_lambda should be 1 in ti0.mdp, -1 in ti1.mdp
        init_lambda should be 0, 1 in ti0.mdp, ti1.mdp
        Set self.current_cycle based on the self.folder_start
        If it's a new start,
            csv file should not exist.
            tpr files should exist.
        If it's a restart,
            csv file should exist. The last line of csv should be consistent with the last cycle (self.folder_start).
            cpt, gro files should exist.
        :return: bool, True if all checks pass
        """
        mdp_list = ["eq0.mdp", "eq1.mdp", "ti0.mdp", "ti1.mdp"]

        # all mdp files should exist
        for mdp in mdp_list:
            if not (self.mdp_folder/mdp).exists():
                logging.error(f"File {self.mdp_folder/mdp} not found")
                return False

        # temperature should be the same
        ref_t_list = [util.get_ref_T(self.mdp_folder / mdp_name) for mdp_name in mdp_list]
        if not np.allclose(ref_t_list, ref_t_list[0], atol=1e-5):
            logging.error(f"Temperature is not the same in all 4 mdp files: {ref_t_list}")
            return False
        self.ref_t = ref_t_list[0]
        self.kBT = util.kB_kj_mol * self.ref_t

        # TI should be 0->1, 1->0
        if not util.mdp_check_TI(self.mdp_folder / "ti0.mdp", 0):
            logging.error(f"nstep * delta_lambda is not 1 in ti0.mdp")
            return False
        if not util.mdp_check_TI(self.mdp_folder / "ti1.mdp", 1):
            logging.error(f"nstep * delta_lambda is not -1 in ti1.mdp")
            return False

        self.current_cycle = int(self.folder_start.name)
        if self.folder_start.name == "000000": # new start
            for tpr in [self.folder_start / '0/eq.tpr', self.folder_start / '1/eq.tpr']:
                if not tpr.exists():
                    logging.error(f"File {tpr} not found")
                    return False
            if self.csv.is_file():
                name_new = util.backup_if_exist_gmx(self.csv)
                logging.info(f"File {self.csv} exists. Backup to {name_new}")
            with open(self.csv, 'w') as f:
                f.writelines([f"Cycle,Work_01 (kJ/mol),Work_10 (kJ/mol),Acceptance_ratio,Accept_{self.ref_t}\n"])
            self.current_folder = self.folder_start

        else: # append run
            if not self.csv.is_file():
                logging.error(f"File {self.csv} not found. Cannot restart.")
                return False
            with open(self.csv) as f:
                lines = f.readlines()
            if len(lines) <= 1:
                logging.error(f"File {self.csv} is empty. Cannot restart")
                return False
            l0 = lines[0].split("_")
            if lines[0][:64] != "Cycle,Work_01 (kJ/mol),Work_10 (kJ/mol),Acceptance_ratio,Accept_":
                logging.error(f"Heading of {self.csv} is incorrect. Cannot restart")
                return False
            if not np.allclose(float(l0[-1]), self.ref_t):
                logging.error(f"Temperature in the heading of {self.csv} is not the same as in mdp files.")
                return False
            last_cycle = int(lines[-1].split(",")[0])
            if last_cycle != self.current_cycle:
                logging.error(f"Last cycle in {self.csv} is {last_cycle}. But the given folder is {self.folder_start.name}")
                return False
            for i in range(2):
                for fname in ["eq.cpt", "eq.gro"]:
                    f = self.folder_start / str(i) / fname
                    if not f.is_file():
                        logging.error(f"{self.csv} suggests restarting from {self.folder_start}, but {f} cannot be found. Restart Failed")
                        return False
            if lines[-1].split(",")[-1] == "A\n":
                self.s0 = [self.folder_start / "1" / "ti.cpt", self.folder_start / "1" / "ti.gro"]
                self.s1 = [self.folder_start / "0" / "ti.cpt", self.folder_start / "0" / "ti.gro"]
                for f in self.s0 + self.s1:
                    if not f.exists():
                        logging.error(f"{self.csv} suggests restarting from {self.folder_start}, but {f} cannot be found. Restart Failed")
                        return False
            elif lines[-1].split(",")[-1] == "R\n":
                self.s0 = [self.folder_start / "0" / "eq.cpt", self.folder_start / "0" / "eq.gro"]
                self.s1 = [self.folder_start / "1" / "eq.cpt", self.folder_start / "1" / "eq.gro"]
            else:
                logging.error(f"Cannot interpreter the last line of {self.csv} : {lines[-1]}")
                return False
            logging.info(
                f"In this restart : 0 starts from {self.s0[0]} {self.s0[1]}, and 1 starts from {self.s1[0]} {self.s1[1]}")
            self.current_cycle += 1
            self.set_current_folder()

        self.safe_flag = True
        self.prepare_scratch_folder()
        if self.min_output:
            if self.folder_start.name == "000000":
                for i in range(2):
                    logging.debug(f"cp {self.folder_start / str(i) / 'eq.tpr'} {self.tmp_folder / str(i) / 'eq.tpr'}")
                    shutil.copy(self.folder_start / str(i) / "eq.tpr", self.tmp_folder / str(i) / "eq.tpr")
        return True

    def log_sim_settings(self):
        """
        Log the simulation settings
        :return: None
        """
        logging.info(f"topology   : {self.top}")
        logging.info(f"csv output : {self.csv}")
        logging.info(f"mdp_folder : {self.mdp_folder}")
        logging.info(f"mdp files  : {self.mdp_folder / 'eq0.mdp'}, {self.mdp_folder / 'eq1.mdp'}, {self.mdp_folder / 'ti0.mdp'}, {self.mdp_folder / 'ti1.mdp'}")
        logging.info(f"Temperature: {self.ref_t} K")
        logging.info(f"kBT        : {self.kBT} kJ/mol")
        logging.info(f"MDRUN      : {self.MDRUN}")
        logging.info(f"GROMPP     : {self.GROMPP}")
        logging.info(f"tmp_folder : {self.tmp_folder}")
        logging.info(f"min_output : {self.min_output}")

    @time_function("run_eq_grompp")
    def run_eq_grompp(self):
        """
        Call self.GROMPP (gmx grompp) to prepare eq.tpr files
        :return: None
        """
        mdp_eq0 = self.mdp_folder / "eq0.mdp"
        mdp_eq1 = self.mdp_folder / "eq1.mdp"

        self.prepare_current_folder()
        cpt0, gro0 = self.s0
        cpt1, gro1 = self.s1
        top_file = self.top
        cmd_list_base = [
            f"{self.GROMPP} -f {mdp_eq0} -c {gro0} -t {cpt0} -p {top_file} ",
            f"{self.GROMPP} -f {mdp_eq1} -c {gro1} -t {cpt1} -p {top_file} ",]
        if self.min_output:
            wdir = self.tmp_folder
            cmd_list = [cmd_list_base[0] + f" -o {wdir / '0/eq.tpr'} -po {wdir/'0/mdout.mdp'} > /dev/null 2>&1",
                        cmd_list_base[1] + f" -o {wdir / '1/eq.tpr'} -po {wdir/'1/mdout.mdp'} > /dev/null 2>&1"]
        else:
            wdir = self.current_folder
            cmd_list = [cmd_list_base[0] + f" -o {wdir / '0/eq.tpr'} -po {wdir/'0/mdout.mdp'} > {wdir / '0' / 'grompp_eq.log'} 2>&1",
                        cmd_list_base[1] + f" -o {wdir / '1/eq.tpr'} -po {wdir/'1/mdout.mdp'} > {wdir / '1' / 'grompp_eq.log'} 2>&1"]
        for cmd in cmd_list:
            logging.debug(cmd)
        for tpr in [wdir / '0/eq.tpr', wdir / '1/eq.tpr']:
            if tpr.exists():
                logging.debug(f"rm {tpr}")
                tpr.unlink()
        processes = [subprocess.Popen(cmd, shell=True) for cmd in cmd_list]
        for p in processes:
            p.communicate()
        # make sure 2 tpr files are generated
        for tpr in [wdir / '0/eq.tpr', wdir / '1/eq.tpr']:
            if not tpr.exists():
                logging.error(f"File {tpr} not found. grompp_eq failed.")
                raise RuntimeError(f"File {tpr} not found. grompp_eq failed.")

    @time_function("run_eq_mdrun")
    def run_eq_mdrun(self):
        """
        Call self.MDRUN (mpirun -np 2 gmx_mpi mdrun) to run eq simulation
        :return: None
        """
        if self.min_output:
            wdir = self.tmp_folder
            cmd = f"{self.MDRUN} -s eq.tpr -multidir {wdir}/0 {wdir}/1 -deffnm eq > /dev/null 2>&1"
        else:
            wdir = self.current_folder
            cmd = f"{self.MDRUN} -s eq.tpr -multidir {wdir}/0 {wdir}/1 -deffnm eq > {wdir / 'mdrun_eq.log'} 2>&1"
        logging.debug(cmd)
        for f in [wdir / "0" / "eq.gro", wdir / "1" / "eq.gro"]:
            if f.exists():
                logging.debug(f"rm {f}")
                f.unlink()
        p = subprocess.Popen(cmd, shell=True, env=self.env)
        p.communicate()
        # make sure 2 eq.gro files are generated
        for f in [wdir / "0" / "eq.gro", wdir / "1" / "eq.gro"]:
            if not f.exists():
                logging.error(f"File {f} not found. eq_mdrun failed.")
                raise RuntimeError(f"File {f} not found. eq_mdrun failed.")
        if self.min_output:
            for f in ["eq.cpt", "eq.gro"]:
                for i in range(2):
                    logging.debug(f"cp {self.tmp_folder / str(i) / f} {self.current_folder / str(i) / f}")
                    shutil.copy(self.tmp_folder / str(i) / f, self.current_folder / str(i) / f)

    @time_function("run_ti_grompp")
    def run_ti_grompp(self):
        """
        Call self.GROMPP (gmx grompp) to prepare ti.tpr files under self.tmp_folder
        :return: None
        """
        mdp_ti0 = self.mdp_folder / "ti0.mdp"
        mdp_ti1 = self.mdp_folder / "ti1.mdp"
        if self.min_output:
            wdir = self.tmp_folder
        else:
            wdir = self.current_folder
        tmp_folder = self.tmp_folder
        cpt0 = wdir / "0" / "eq.cpt"
        gro0 = wdir / "0" / "eq.gro"
        tpr0 = tmp_folder / "0" / "ti.tpr"
        mdout0 = tmp_folder / "0" / "mdout.mdp"
        cpt1 = wdir / "1" / "eq.cpt"
        gro1 = wdir / "1" / "eq.gro"
        tpr1 = tmp_folder / "1" / "ti.tpr"
        mdout1 = tmp_folder / "1" / "mdout.mdp"
        top_file = self.top
        cmd_list_base = [
            f"{self.GROMPP} -f {mdp_ti0} -c {gro0} -t {cpt0} -p {top_file} -o {tpr0} -po {mdout0} ",
            f"{self.GROMPP} -f {mdp_ti1} -c {gro1} -t {cpt1} -p {top_file} -o {tpr1} -po {mdout1} ",
        ]
        if self.min_output:
            cmd_list = [cmd + " > /dev/null 2>&1" for cmd in cmd_list_base]
        else:
            cmd_list = [cmd_list_base[0] + f" > {wdir / '0' / 'grompp_ti.log'} 2>&1 ",
                        cmd_list_base[1] + f" > {wdir / '1' / 'grompp_ti.log'} 2>&1 ",
                        ]
        for cmd in cmd_list:
            logging.debug(cmd)
        for tpr in [tpr0, tpr1]:
            if tpr.exists():
                logging.debug(f"rm {tpr}")
                tpr.unlink()
        processes = [subprocess.Popen(cmd, shell=True) for cmd in cmd_list]
        for p in processes:
            p.communicate()
        # make sure the tpr files are generated
        for tpr in [tpr0, tpr1]:
            if not tpr.exists():
                logging.error(f"File {tpr} not found. grompp_ti failed.")
                raise RuntimeError(f"File {tpr} not found. grompp_ti failed.")

    @time_function("run_ti_mdrun")
    def run_ti_mdrun(self):
        """
        Call self.MDRUN (mpirun -np 2 gmx_mpi mdrun) to run TI simulation
        :return: None
        """
        mdrun = self.MDRUN
        tmp_folder = self.tmp_folder
        wdir = self.current_folder
        multi_dir = f"-multidir {tmp_folder}/0 {tmp_folder}/1"
        cmd = f"{mdrun} -s ti.tpr {multi_dir} -deffnm ti "
        if self.min_output:
            cmd += " > /dev/null 2>&1"
        else:
            cmd += f" > {wdir / 'mdrun_ti.log'} 2>&1"
        logging.debug(cmd)
        for f in [tmp_folder / "0" / "ti.gro", tmp_folder / "1" / "ti.gro"]:
            if f.exists():
                logging.debug(f"rm {f}")
                f.unlink()
        p = subprocess.Popen(cmd, shell=True, env=self.env)
        p.communicate()
        # make sure 2 ti.gro files are generated
        for f in [tmp_folder / "0" / "ti.gro", tmp_folder / "1" / "ti.gro"]:
            if not f.exists():
                logging.error(f"File {f} not found. ti_mdrun failed.")
                raise RuntimeError(f"File {f} not found. ti_mdrun failed.")

        # integrate. get work value. W<0 means system releases energy
        work01 = util.integrate_work(tmp_folder / "0" / "ti.xvg")
        work10 = -util.integrate_work(tmp_folder / "1" / "ti.xvg")
        swap_flag, csv_line = self.swap_check(work01, work10)
        if swap_flag: # save ti.gro and ti.cpt
            for i in range(2):
                tmp_folder = self.tmp_folder / str(i)
                rep_folder = self.current_folder / str(i)
                for f in ["ti.gro", "ti.cpt"]:
                    shutil.copy(tmp_folder / f, rep_folder / f)

        if self.debug: # if debug model, save ti.xvg ti.tpr
            for i in range(2):
                tmp_folder = self.tmp_folder / str(i)
                rep_folder = self.current_folder / str(i)
                for f in ["ti.xvg", "ti.tpr"]:
                    shutil.copy(tmp_folder / f, rep_folder / f)
        # clean up, except ti.gro, ti.cpt
        for f in [self.tmp_folder/"0", self.tmp_folder/"1"]:
            for f2 in f.iterdir():
                if f2.is_file():
                    if f2.name not in ["ti.gro", "ti.cpt", "eq.cpt", "eq.gro"]:
                        logging.debug(f"rm {f2}")
                        f2.unlink()
        
        with open(self.csv, 'a') as f:
            f.write(csv_line)
        

    def swap_check(self, w01, w10):
        """

        :param w01: Work 0->1 in kJ/mol
        :param w10: Work 1->0 in kJ/mol
        :return: swap_flag, csv_line
        """
        p_accept = np.exp(-(w01 + w10) / self.kBT)
        csv_line = f"{self.current_cycle:5d},{w01:16.12},{w10:16.12},{min(1, p_accept):16.3f},"
        info_line = f"Cycle {self.current_cycle}, {w01:9.3f} kJ/mol, {w10:9.3f} kJ/mol, {min(1, p_accept):6.3f},"
        swap_flag = np.random.rand() < p_accept
        if swap_flag:
            csv_line += "A\n"
            info_line += " Accept"
            self.s0 = [self.tmp_folder / "1" / "ti.cpt", self.tmp_folder / "1" / "ti.gro"]
            self.s1 = [self.tmp_folder / "0" / "ti.cpt", self.tmp_folder / "0" / "ti.gro"]
        else:
            csv_line += "R\n"
            info_line += " Reject"
            if self.min_output:
                self.s0 = [self.tmp_folder / "0" / "eq.cpt", self.tmp_folder / "0" / "eq.gro"]
                self.s1 = [self.tmp_folder / "1" / "eq.cpt", self.tmp_folder / "1" / "eq.gro"]
            else:
                self.s0 = [self.current_folder / "0" / "eq.cpt", self.current_folder / "0" / "eq.gro"]
                self.s1 = [self.current_folder / "1" / "eq.cpt", self.current_folder / "1" / "eq.gro"]
        logging.info(info_line)
        logging.info(f"In the next cycle: 0 starts from {self.s0[0]} {self.s0[1]}, and 1 starts from {self.s1[0]} {self.s1[1]}")

        # with open(self.csv, 'a') as f:
        #     f.write(csv_line)
        return swap_flag, csv_line

    def cycle(self, cycle, maxh=24):
        """
        :param cycle: int, number of cycles (work evaluation) to run
        :param maxh: float, terminate after this time. Time will only be checked at the start of a cycle.
        If either cycle or maxh is reached, the program will stop.
        :return: Bool, True if the simulation is finished without error.
        """
        if not self.safe_flag:
            logging.error("Safety check failed. Cannot start the simulation.")
            return False
        t0 = time.time()
        for i in range(cycle):
            t1 = time.time()
            if t1 - t0 > maxh * 3600:
                logging.info(f"Terminate after {maxh} h. Time exceeded.")
                return False
            logging.info(f"Cycle {self.current_cycle}, eq")
            if self.current_cycle != 0:
                self.run_eq_grompp()
            self.run_eq_mdrun()
            logging.info(f"Cycle {self.current_cycle}, TI")
            self.run_ti_grompp()
            self.run_ti_mdrun()
            self.current_cycle += 1
            self.set_current_folder()

            t2 = time.time()
            logging_ave_time_cycle((t2-t0)/(i+1) , t2 - t1)
        return True

    def estimate_free_energy(self):
        """
        Read csv file and estimate free energy using BAR
        :return: dG, dG_error
        """
        df = pd.read_csv(self.csv)
        work01 = df["Work_01 (kJ/mol)"].values
        work10 = df["Work_10 (kJ/mol)"].values
        dG, dGe = util.free_E_bar(work01/self.kBT, work10/self.kBT)
        return dG, dGe, df


def main():
    t0 = time.time()
    parser = argparse.ArgumentParser(
        description=f"""Version {pmxNCMC.__version__}. 
            This is the IO based implementation for running pmx in expanded ensemble. It can also be understand as 
            Replica Exchange with lambda 0 and 1 only or Non-equilibrium Candidate Monte Carlo.""")
    parser.add_argument('-p',
                        metavar='topology',
                        type=Path, help='Topology file',
                        default='topol.top')
    parser.add_argument('-log', metavar=" ",
                        type=str, help='Log file',
                        default='md.log')
    parser.add_argument('-csv', metavar=" ",
                        type=Path, help='CSV file for saving work values. Default : md.csv',
                        default='md.csv')
    parser.add_argument('-mdp_folder', metavar=" ",
                        type=Path, default=Path("mdp"),
                        help='Folder that contains eq0.mdp, eq1.mdp, ti0.mdp, ti1.mdp. All 4 files are required. '
                             'File name should be exact. Default : ./mdp')
    parser.add_argument('-folder_start', metavar=" ",
                        type=str, help='Folder to start the simulation. If not given, csv file will be checked. ')
    parser.add_argument('-cycle', metavar=" ",
                        type=lambda x: int(x) if int(x) > 0 else parser.error("cycle must be greater than 0"),
                        help='Number of cycles (work evaluation) to run.')
    parser.add_argument('-cyc_until', metavar=" ",
                        type=lambda x: int(x) if int(x) > 0 else parser.error("cyc_until must be greater than 0"),
                        help='Number of cycles (work evaluation) to run.',)
    parser.add_argument('-maxh', metavar=" ",
                        type=float, help='Terminate after this time. Time will only be checked at the start of a cycle. '
                                         'The actually running time can possibly exceed this time. Default : 23.5 h',
                        default=23.5)
    parser.add_argument('-MDRUN', metavar=" ",
                        type=str,
                        help='Command for mdrun, we will use multidir, MPI is required. '
                             'Default : "mpirun -np 2 --bind-to none gmx_mpi mdrun"',
                        default='mpirun -np 2 --bind-to none gmx_mpi mdrun')
    parser.add_argument('-GROMPP', metavar=" ",
                        type=str, help='Command for grompp, with additional flags. '
                                       'For example "gmx_threads_AVX_256 grompp -maxwarn 1" Default : "gmx grompp"',
                        default='gmx grompp')
    parser.add_argument('-tmp_folder', metavar=" ",
                        type=str,
                        help='Temporary folder. Point it to the local storage on the computing node to save IO. '
                             'You can also set this to /dev/shm if you have enough memory. '
                             'Default : auto determined by python tempfile. https://docs.python.org/3.11/library/tempfile.html', )
    parser.add_argument('--debug',
                        action='store_true', help='Print debug information and save all TI output.')
    parser.add_argument('--format', metavar=" ",
                        type=str,
                        help='Log format. Default : "%%(message)s" . In debug run, '
                             'the format will be forced to "%%(asctime)s - %%(levelname)s - %%(message)s" . '
                             'Please check https://docs.python.org/3/library/logging.html#logging.basicConfig for '
                             'more detail.',
                        default='%(message)s')
    parser.add_argument('--min_output',
                        action='store_true',
                        help='Minimal IO. If you are sure about what you are doing and you want the minimal output '
                             'to save IO. This will run eq on tmp_folder and redirect all gromacs stdout/stderr to '
                             '/dev/null. Only eq.cpt, eq.gro, ti.cpt, ti.gro will be saved.')

    args = parser.parse_args()
    env = os.environ.copy()

    if not (args.folder_start is None):
        folder_start = args.folder_start
    else:
        if args.csv.is_file():
            with open(args.csv) as f:
                lines = f.readlines()
            if len(lines) > 1:
                current_cycle = int(lines[-1].split(",")[0])
                folder_start = f"{current_cycle:06d}"
            else:
                folder_start = "000000"
        else:
            folder_start = "000000"
    if folder_start == "000000":  # new start
        util.backup_if_exist_gmx(args.log)

    if args.debug:
        logging.basicConfig(
            filename=args.log, filemode='a',
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s')
        logging.debug(f"ti.xvg and ti.tpr will be saved.")
    else:
        logging.basicConfig(
            filename=args.log, filemode='a',
            level=logging.INFO,
            format=args.format)
    command_line = ""

    for word in sys.argv:
        if " " in word or "%" in word:
            command_line += f'"{word}" '
        else:
            command_line += word + " "
    logging.info(f"Command line: {command_line}")
    logging.info(f"pmx_mdrun version {pmxNCMC.__version__}")
    logging.info(f"log file: {args.log}")

    if args.tmp_folder is None:
        tmp_folder = Path(tempfile.mkdtemp(prefix="pmxRE_"))
    else:
        tmp_folder = Path(tempfile.mkdtemp(dir=args.tmp_folder, prefix="pmxRE_"))



    # count cycle if cyc_until is given
    if args.cyc_until is not None:
        if args.cycle is not None:
            logging.info("Argument cycle will be ignored.")
        args.cycle = args.cyc_until - int(folder_start) - 1
        if int(folder_start) == 0:
            args.cycle += 1
        if args.cycle <= 0:
            logging.error(f"Cannot start from {folder_start} and run until {args.cyc_until}")
            return
        else:
            logging.info(f"Run until cycle {args.cyc_until}. {args.cycle} cycles to run.")
    else:
        if args.cycle is None:
            parser.error("Either cycle or cyc_until must be given.")



    try:
        mdrun = PMX_MDRUN_RE(top=args.p,
                             csv=args.csv,
                             mdp_folder=args.mdp_folder,
                             folder_start=Path(folder_start),
                             MDRUN=args.MDRUN,
                             GROMPP=args.GROMPP,
                             tmp_folder=tmp_folder,
                             env=env,
                             min_output=args.min_output,
                             debug=args.debug)
        if not mdrun.safety_check():
            return
        mdrun.log_sim_settings()
        s_flag = mdrun.cycle(args.cycle, args.maxh)
        if s_flag:
            dG, dGe, df = mdrun.estimate_free_energy()
            logging.info(f"Estimate free energy difference using {len(df)} cycles")
            kBT = mdrun.kBT
            logging.info(f"DeltaG = {dG * kBT:6.2f} +- {dGe * kBT:4.2f} kJ/mol")
            kBT_kcal = util.kB_kcal_mol * mdrun.ref_t
            logging.info(f"       = {dG * kBT_kcal:6.2f} +- {dGe * kBT_kcal:4.2f} kcal/mol")

            # log the time for each function
            t1 = time.time() - t0
            for key, value in fun_exe_times.items():
                if len(value) > 0:
                    logging.info(f"Average time {key:13s} :{np.mean(value):7.2f} s, {np.sum(value)/t1*100:5.2f}%")
    finally:
        shutil.rmtree(tmp_folder)
        t1 = time.time() - t0
        logging.info(f"pmx_mdrun finished in {int(t1 // 3600)} h {int((t1 % 3600) // 60)} min {(t1 % 60):.1f} s")



if __name__ == "__main__":
    main()

