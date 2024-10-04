#!/usr/bin/env python3
import argparse
import sys
import shutil
from pathlib import Path
import time
import logging
import tempfile

import numpy as np
from scipy import integrate

import pmxNCMC
from pmxNCMC import util
import subprocess
from concurrent.futures import ThreadPoolExecutor


def prepare_scratch_folder(tmp_folder_path):
    if not tmp_folder_path.exists():
        tmp_folder_path.mkdir()
    for i in range(2):
        tmp_ = tmp_folder_path / str(i)
        if not tmp_.exists():
            tmp_.mkdir()
    # make sure the folder is empty
    for i in range(2):
        tmp_ = tmp_folder_path / str(i)
        for f in tmp_.iterdir():
            if f.is_file():
                f.unlink()
            else:
                shutil.rmtree(f)

def prepare_current_folder(current_folder_path):
    for i in range(2):
        rep_dir = current_folder_path / str(i)
        if not rep_dir.exists():
            rep_dir.mkdir(parents=True)

def run_eq_grompp(settings_dict, s0, s1):
    """
    Run eq simulation
    gmx grompp
    mpirun -np 2 gmx_mpi mdrun -multidir
    """
    mdp_eq0 = settings_dict["mdp_folder"] / "eq0.mdp"
    mdp_eq1 = settings_dict["mdp_folder"] / "eq1.mdp"
    current_folder = settings_dict["current_folder"]
    cpt0, gro0 = s0
    cpt1, gro1 = s1
    top_file = settings_dict["top"]
    command_list = [
        f"{settings_dict['GROMPP']} -f {mdp_eq0} -c {gro0} -t {cpt0} -p {top_file} -o {current_folder / '0' / 'eq.tpr'} > {current_folder / '0' / 'grompp_eq.log'} 2>&1",
        f"{settings_dict['GROMPP']} -f {mdp_eq1} -c {gro1} -t {cpt1} -p {top_file} -o {current_folder / '1' / 'eq.tpr'} > {current_folder / '1' / 'grompp_eq.log'} 2>&1",
    ]
    for cmd in command_list:
        logging.debug(cmd)
    processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for cmd in command_list]

    # Collect the results as the commands finish
    for p in processes:
        p.wait()
    # make sure the tpr files are generated
    for tpr in [current_folder / '0' / 'eq.tpr', current_folder / '1' / 'eq.tpr']:
        if not tpr.exists():
            logging.info(f"File {tpr} not found")
            raise RuntimeError(f"File {tpr} not found")

def run_eq_mdrun(settings_dict):
    """
    Run eq simulation.
    Assume all tpr file has been generated. current_folder/0/eq.tpr, current_folder/1/eq.tpr
    :param settings_dict:
    :return: None
    """
    current_folder = settings_dict["current_folder"]
    mdrun = settings_dict["MDRUN"]
    multi_dir = f"-multidir {current_folder}/0 {current_folder}/1"
    command = f"{mdrun} -s eq.tpr {multi_dir} -deffnm eq > {current_folder / 'mdrun_eq.log'} 2>&1"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    # make sure 2 eq.gro files are generated
    for f in [current_folder / "0" / "eq.gro", current_folder / "1" / "eq.gro"]:
        if not f.exists():
            raise RuntimeError(f"After eq_mdrun, {f} was not found.")

def run_ti_grompp(settings_dict):
    """
    grompp for TI simulation, 0->1 starts from 0 (cpt, gro), 1->0 starts from 1 (cpt, gro)
    :return: None
    """
    mdp_ti0 = settings_dict["mdp_folder"] / "ti0.mdp"
    mdp_ti1 = settings_dict["mdp_folder"] / "ti1.mdp"
    current_folder = settings_dict["current_folder"]
    tmp_folder = settings_dict["tmp_folder"]
    cpt0 = current_folder / "0" / "eq.cpt"
    gro0 = current_folder / "0" / "eq.gro"
    tpr0 = tmp_folder / "0" / "ti.tpr"
    cpt1 = current_folder / "1" / "eq.cpt"
    gro1 = current_folder / "1" / "eq.gro"
    tpr1 = tmp_folder / "1" / "ti.tpr"
    top_file = settings_dict["top"]
    command_list = [
        f"{settings_dict['GROMPP']} -f {mdp_ti0} -c {gro0} -t {cpt0} -p {top_file} -o {tpr0} > {current_folder / '0' / 'grompp_ti.log'} 2>&1",
        f"{settings_dict['GROMPP']} -f {mdp_ti1} -c {gro1} -t {cpt1} -p {top_file} -o {tpr1} > {current_folder / '1' / 'grompp_ti.log'} 2>&1",
    ]
    for cmd in command_list:
        logging.debug(cmd)
    processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for cmd in command_list]
    # Collect the results as the commands finish
    for p in processes:
        p.wait()
    # make sure the tpr files are generated
    for tpr in [tpr0, tpr1]:
        if not tpr.exists():
            raise RuntimeError(f"After ti_grompp {tpr} was not found.")

def run_ti_mdrun(settings_dict):
    """
    Run TI simulation, return work value in kJ/mol
    :param settings_dict:
    :return: work01, work10
    """
    mdrun = settings_dict["MDRUN"]
    tmp_folder = settings_dict["tmp_folder"]
    current_folder = settings_dict["current_folder"]
    multi_dir = f"-multidir {tmp_folder}/0 {tmp_folder}/1"
    command = f"{mdrun} -s ti.tpr {multi_dir} -deffnm ti > {current_folder / 'mdrun_ti.log'} 2>&1"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    # make sure 2 ti.gro files are generated
    for f in [tmp_folder / "0" / "ti.gro", tmp_folder / "1" / "ti.gro"]:
        if not f.exists():
            raise RuntimeError(f"After ti_mdrun, {f} not found.")

    # integrate. get work value. The work that external force does on the system. <0 means system releases energy
    work01 =  util.integrate_work(tmp_folder / "0" / "ti.xvg")
    work10 = -util.integrate_work(tmp_folder / "1" / "ti.xvg")

    # save ti.gro and ti.cpt
    for i in range(2):
        tmp_folder = settings_dict["tmp_folder"] / str(i)
        rep_folder = settings_dict["current_folder"] / str(i)
        for f in ["ti.gro", "ti.cpt"]:
            shutil.copy(tmp_folder / f, rep_folder / f)
    # if dry_run, also save ti.xvg ti.tpr
    if settings_dict["DEBUG"]:
        for i in range(2):
            tmp_folder = settings_dict["tmp_folder"] / str(i)
            rep_folder = settings_dict["current_folder"] / str(i)
            for f in ["ti.xvg", "ti.tpr"]:
                shutil.copy(tmp_folder / f, rep_folder / f)
    # clean up, except ti.gro, ti.cpt
    for f in [settings_dict["tmp_folder"]/"0", settings_dict["tmp_folder"]/"1"]:
        for f2 in f.iterdir():
            if f2.is_file():
                if f2.name not in ["ti.gro", "ti.cpt"]:
                    logging.debug(f"Remove file {f2}")
                    f2.unlink()
            else:
                shutil.rmtree(f2)
                logging.debug(f"Remove folder {f2}")
    return work01, work10

def swap_check(w01, w10, kBT, settings_dict):
    """
    :param w01: work 0->1 in kJ/mol
    :param w10: work 1->0 in kJ/mol
    :param kBT: kJ/mol
    :param settings_dict:
    Work is defined as the work that external force does on the system. <0 means system releases energy
    :return: swap_flag, csv_line, s0, s1
    """
    p_accept = np.exp(-(w01 + w10) / kBT)
    csv_line = f"{settings_dict['current_cycle']:5d},{w01:16.12},{w10:16.12},{min(1, p_accept):16.3f},"
    inf_line = f"Cycle {settings_dict['current_cycle']}, {w01:9.3f} kJ/mol, {w10:9.3f} kJ/mol, {min(1, p_accept):6.3f},"
    current_folder = settings_dict["current_folder"]
    tmp_folder = settings_dict["tmp_folder"]
    swap_flag = np.random.rand() < p_accept
    if swap_flag:
        csv_line += "A\n"
        inf_line += " Accept"
        s0 = tmp_folder / "1/ti.cpt", tmp_folder / "1/ti.gro"
        s1 = tmp_folder / "0/ti.cpt", tmp_folder / "0/ti.gro"
    else:
        csv_line += "R\n"
        inf_line += " Reject"
        s0 = current_folder / "0/eq.cpt", current_folder / "0/eq.gro"
        s1 = current_folder / "1/eq.cpt", current_folder / "1/eq.gro"
    logging.info(inf_line)
    logging.debug(f"In the next cycle: 0 starts from {s0[0]} {s0[1]}, and 1 starts from {s1[0]} {s1[1]}")
    return swap_flag, csv_line, s0, s1

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


if __name__ == "__main__":
    t0 = time.time()
    parser = argparse.ArgumentParser(
        description = f"""Version {pmxNCMC.__version__}. 
        This is the IO based implementation for running pmx in expanded ensemble. It can also be understand as 
        Replica Exchange with lambda 0 and 1 only or Non-equilibrium Candidate Monte Carlo.""")
    parser.add_argument('-p',
                        metavar='topology',
                        type=str, help='Topology file',
                        default='topol.top')
    parser.add_argument('-log',
                        type=str, help='Log file',
                        default='md.log')
    parser.add_argument('-csv',
                        metavar='csv output',
                        type=str, help='CSV file for saving work values. Default : md.csv',
                        default='md.csv')
    parser.add_argument('-mdp_folder',
                        type=str,
                        help='Folder that contains eq0.mdp, eq1.mdp, ti0.mdp, ti1.mdp. All 4 files are required. '
                             'File name should be exact.')
    parser.add_argument('-folder_start',
                        type=str, help='Folder to start the simulation. Default : 000000',
                        default='000000')
    parser.add_argument('-cycle',
                        type=int, help='Number of cycles (work evaluation) to run. Default : 10',
                        default=10)
    parser.add_argument('-maxh',
                        type=float, help='Terminate after this time. It will only be checked at the start of a cycle. '
                                         'The actually running time can possibly exceed this time. Default : 23.5 h',
                        default=23.5)
    parser.add_argument('-MDRUN', metavar="",
                        type=str, help='command for mdrun, we will use multidir, MPI is required. Default : "mpirun -np 2 gmx_mpi mdrun"',
                        default='mpirun -np 2 --bind-to none gmx_mpi mdrun')
    parser.add_argument('-GROMPP', metavar="",
                        type=str, help='command for grompp, with additional flags. '
                                       'For example "gmx_threads_AVX_256 grompp -maxwarn 1" Default : "gmx grompp"',
                        default='gmx grompp')
    parser.add_argument('-tmp_folder',
                        type=str, help='Temporary folder. Point it to the local storage on the computing node to save IO. '
                                       'You can also set this to /dev/shm if you have enough memory. '
                                       'Default : auto determined by python tempfile',)
    parser.add_argument('-re_try', help='Number of re-try if the simulation fails. Default : 3',
                        type=int,
                        default=3)
    parser.add_argument('--debug',
                        action='store_true', help='Print debug information')
    parser.add_argument('--format', metavar="log_format",
                        type=str, help='Log format. Default : "%%(asctime)s - %%(levelname)s - %%(message)s" . '
                                       'If you want a clean output, use "%%(message)s". In debug run, this will be ignored.',
                        default='%(asctime)s - %(levelname)s - %(message)s')

    args = parser.parse_args()
    settings = {"top"          : Path(args.p),
                "log"          : Path(args.log),
                "csv"          : Path(args.csv),
                "mdp_folder"   : Path(args.mdp_folder),
                "cycle"        : args.cycle,
                "folder_start" : Path(args.folder_start),
                "maxh"         : args.maxh,
                "MDRUN"        : args.MDRUN,
                "GROMPP"       : args.GROMPP,
                "tmp_folder"   : args.tmp_folder,
                "DEBUG"        : args.debug,
                "current_cycle": 0,
                "base_path"    : Path.cwd(),
                "re_try"       : args.re_try,
                }
    if settings["log"].exists():
        util.backup_if_exist(settings["log"])
    if args.debug:
        # set log to DEBUG
        logging.basicConfig(
            filename=args.log, filemode='w',
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s')
        logging.debug(f"ti.xvg and ti.tpr will be saved.")
    else:
        # print(args.log, args.format)
        logging.basicConfig(
            filename=args.log, filemode='w',
            level=logging.INFO,
            format=args.format)


    if args.tmp_folder is None:
        settings["tmp_folder"] = Path(tempfile.mkdtemp(prefix="pmxNCMCRE_"))
    else:
        settings["tmp_folder"] = Path(tempfile.mkdtemp(dir=args.tmp_folder, prefix="pmxNCMCRE_"))

    # make sure input files exist
    for name in ["top"]:
        if not settings[name].exists():
            logging.info(f"File {settings[name]} not found")
            exit(1)

    for mdp_name in ["eq0.mdp", "eq1.mdp", "ti0.mdp", "ti1.mdp"]:
        if not (settings["mdp_folder"] / mdp_name).exists():
            logging.info(f"File {settings['mdp_folder'] / mdp_name} not found")
            exit(1)
    ref_t_list = [util.get_ref_T(settings["mdp_folder"] / mdp_name) for mdp_name in ["eq0.mdp", "eq1.mdp", "ti0.mdp", "ti1.mdp"]]
    # all close
    if not np.allclose(ref_t_list, ref_t_list[0]):
        logging.info(f"Reference temperature are not the same: {ref_t_list}")
        exit(1)
    else:
        settings["ref_t"] = ref_t_list[0]
    for mdp in ["ti0.mdp", "ti1.mdp"]:
        if not util.mdp_check_TI(settings["mdp_folder"] / mdp):
            exit(1)

    # if folder_start is 00000, we will start from scratch
    if settings["folder_start"].name == "000000":
        logging.info("New simulation starting from cycle 0")
        settings["current_cycle"] = 0
        for file_name in [settings["csv"]]:
            if file_name.exists():
                util.backup_if_exist(file_name)
        with open(settings["csv"], "w") as f:
            f.writelines([f"Cycle,Work_01 (kJ/mol),Work_10 (kJ/mol),Acceptance_ratio,Accept_{settings['ref_t']}\n"])
    elif int(settings["folder_start"].name) < 0:
        logging.info(f"Invalid folder_start {settings['folder_start']}")
        exit(1)
    else:
        settings["current_cycle"] = int(settings["folder_start"].name)
        logging.info(f"Restart from Cycle {settings['current_cycle']}, assume this Cycle has been finished")
        # read csv, prepare s0, s1
        with open(settings["csv"], "r") as f:
            l = f.readlines()[-1].rstrip()
        words = l.split(",")
        cycle = int(words[0])
        # w01, w10, p_accept= [float(i) for i in words[1:-1]]
        # check if cycle from csv matches the folder
        if cycle != settings["current_cycle"]:
            logging.info(f"Cycle from csv {cycle} does not match the folder {settings['current_cycle']}")
            exit(1)
        else:
            current_folder = Path(f"{settings['current_cycle']:06d}")

        # prepare starting cpt,gro
        if words[-1] == "A":
            s0 = current_folder / "1/ti.cpt", current_folder / "1/ti.gro"
            s1 = current_folder / "0/ti.cpt", current_folder / "0/ti.gro"
        elif words[-1] == "R":
            s0 = current_folder / "0/eq.cpt", current_folder / "0/eq.gro"
            s1 = current_folder / "1/eq.cpt", current_folder / "1/eq.gro"
        else:
            logging.info(f"Invalid csv line {l}, cannot proceed with restart.")
            exit(1)
        for f in s0 + s1:
            logging.debug(f"Check file {f}")
            if not f.exists():
                logging.info(f"File {f} not found")
                exit(1)
        settings['current_cycle'] += 1
        current_folder = Path(f"{settings['current_cycle']:06d}")

    if settings["re_try"] <= 0:
        logging.info(f"re_try should be larger than 0")
        exit(1)



    kBT = 8.314462618e-3 * settings["ref_t"] # kJ * K/mol
    command_line = ""
    for word in sys.argv:
        if " " in word:
            command_line += f'"{word}" '
        else:
            command_line += word + " "
    logging.info(f"Command line: {command_line}")
    logging.info(f"pmx_mdrun version {pmxNCMC.__version__}")
    logging.info(f"# Simulation settings #############################################################################")
    logging.info(f"topology   : {settings['top']}")
    logging.info(f"log output : {settings['log']}")
    logging.info(f"csv output : {settings['csv']}")
    logging.info(f"mdp_folder : {settings['mdp_folder']}")
    logging.info(f"mdp files  : {settings['mdp_folder']/'eq0.mdp'} {settings['mdp_folder']/'eq1.mdp'} {settings['mdp_folder']/'ti0.mdp'} {settings['mdp_folder']/'ti1.mdp'}")
    logging.info(f"Temperature: {settings['ref_t']} K")
    logging.info(f"kBT        : {kBT} kJ/mol")
    logging.info(f"MDRUN  command : {settings['MDRUN']}")
    logging.info(f"GROMPP command : {settings['GROMPP']}")
    logging.info(f"tmp_folder     : {settings['tmp_folder']}")
    logging.info(f"Cycle (eq+ti) to run : {settings['cycle']}")
    logging.info(f"Maximum running time : {settings['maxh']} h = {settings['maxh']*60} min")
    logging.info(f"Re-try if failed     : {settings['re_try']}")
    logging.info(f"# Simulation settings End #########################################################################")

    prepare_scratch_folder(settings["tmp_folder"])
    try:
        # eq+NCMC cycle
        t_tick = time.time()  # track time for each cycle
        for cycle in range(settings["cycle"]):
            # check time
            if (time.time() - t0) / 3600 > settings["maxh"]:
                logging.info(f"Time limit reached, exit")
                break

            settings["current_folder"] = Path(f"{settings['current_cycle']:06d}")
            logging.info(f"Cycle {settings['current_cycle']}, eq")
            prepare_current_folder(settings['current_folder']) # mkdir current_folder/0, current_folder/1


            succ_flag = False
            for re_try in range(settings["re_try"]):
                try:
                    # eq
                    if settings["current_cycle"] != 0:
                        run_eq_grompp(settings, s0, s1)
                    run_eq_mdrun(settings)
                    # TI
                    logging.info(f"Cycle {settings['current_cycle']}, TI")
                    run_ti_grompp(settings)
                    w01, w10 = run_ti_mdrun(settings)
                    succ_flag = True
                    break
                except RuntimeError as e:
                    logging.info(f"Cycle {settings['current_cycle']} failed, retry {re_try+1}")
                    logging.info(f"Error: {e}")
            if not succ_flag:
                logging.info(f"Cycle {settings['current_cycle']} failed {settings['re_try']} fimes, exit")
                raise RuntimeError(f"Cycle {settings['current_cycle']} failed {settings['re_try']} fimes, exit")

            # swap attempt
            _, csv_line, s0, s1 =  swap_check(w01, w10, kBT, settings)
            with open(settings["csv"], "a") as f:
                f.write(csv_line)

            settings["current_cycle"] += 1
            t_tock = time.time()
            logging_ave_time_cycle((t_tock-t0) / (cycle+1), t_tock-t_tick)
            t_tick = t_tock

    finally:
        t1 = time.time()-t0
        logging.info(f"pmx_mdrun finished in {int(t1 // 3600)} h {int((t1 % 3600) // 60)} min {(t1 % 60):.1f} s")
        logging.info(f"Cleaning up ...")
        shutil.rmtree(settings["tmp_folder"])
        logging.info(f"ALL Done")


