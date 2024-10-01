#!/usr/bin/env python3
import argparse
import os
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
    processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for cmd in command_list]

    # Collect the results as the commands finish
    for p in processes:
        stdout, stderr = p.communicate()
    # make sure the tpr files are generated
    for tpr in [current_folder / '0' / 'eq.tpr', current_folder / '1' / 'eq.tpr']:
        if not tpr.exists():
            logging.info(f"File {tpr} not found")
            raise Exception(f"File {tpr} not found")

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
    processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for cmd in command_list]
    # Collect the results as the commands finish
    for p in processes:
        stdout, stderr = p.communicate()
    # make sure the tpr files are generated
    for tpr in [tpr0, tpr1]:
        if not tpr.exists():
            logging.info(f"File {tpr} not found")
            raise Exception(f"File {tpr} not found")

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

    # get work value, read ti.xvg and integrate, column 1 is time, column 2 is dH/dlambda
    ti01 = np.loadtxt(tmp_folder / "0" / "ti.xvg", comments=["#", "@"])
    ti10 = np.loadtxt(tmp_folder / "1" / "ti.xvg", comments=["#", "@"])
    # integrate
    lam = np.linspace(0, 1, len(ti01))
    work01 = -integrate.simpson(ti01[:, 1], x=lam)
    lam = np.linspace(1, 0, len(ti10))
    work10 = -integrate.simpson(ti10[:, 1], x=lam)

    # save ti.gro and ti.cpt
    for i in range(2):
        tmp_folder = settings_dict["tmp_folder"] / str(i)
        rep_folder = settings_dict["current_folder"] / str(i)
        for f in ["ti.gro", "ti.cpt"]:
            shutil.copy(tmp_folder / f, rep_folder / f)
    # clean up everything in under 0 and 1
    for f in [settings_dict["tmp_folder"]/"0", settings_dict["tmp_folder"]/"1"]:
        for f2 in f.iterdir():
            if f2.is_file():
                f2.unlink()
            else:
                shutil.rmtree(f2)
    return work01, work10

def swap_check(w01, w10, kBT, settings):
    """
    :param w01: work 0->1 in kJ/mol
    :param w10: work 1->0 in kJ/mol
    :param kBT: kJ*K/mol
    :param settings:
    :return: s0, s1
    """
    p_accept = np.exp(-(w01 + w10) / kBT)
    csv_line = f"{settings['current_cycle']},{w01},{w10},{min(1,p_accept):6.3f},"
    inf_line = f"Cycle {settings['current_cycle']},{w01:10.3f},{w10:10.3f}, {min(1,p_accept):6.3f},"
    current_folder = settings["current_folder"]
    if np.random.rand() < p_accept:
        csv_line += "A\n"
        inf_line += " Accept"
        s0 = current_folder / "1/ti.cpt", current_folder / "1/ti.gro"
        s1 = current_folder / "0/ti.cpt", current_folder / "0/ti.gro"
    else:
        csv_line += "R\n"
        inf_line += " Reject"
        s0 = current_folder / "0/ti.cpt", current_folder / "0/ti.gro"
        s1 = current_folder / "1/ti.cpt", current_folder / "1/ti.gro"
    logging.info(inf_line)
    return csv_line, s0, s1


if __name__ == "__main__":
    t0 = time.time()
    parser = argparse.ArgumentParser(
        description = f"""Version {pmxNCMC.__version__}. 
        This is the IO based implementation for running pmx in expanded ensemble. It can also be understand as 
        Replica Exchange with lambda 0 and 1 only or Non-equilibrium Candidate Monte Carlo.""")
    parser.add_argument('-p',
                        type=str, help='Topology file',
                        default='topol.top')
    parser.add_argument('-log',
                        type=str, help='Log file',
                        default='md.log')
    parser.add_argument('-csv',
                        type=str, help='CSV file for saving work values. Default : md.csv',
                        default='md.csv')

    parser.add_argument('-mdp_folder',
                        type=str,
                        help='Folder that contains eq0.mdp, eq1.mdp, ti0.mdp, ti1.mdp. All 4 files are required. Please don\'t save trajectory.')
    parser.add_argument('-folder_start',
                        type=str, help='Folder to start the simulation.', required=True)
    parser.add_argument('-cycle',
                        type=int, help='Number of cycles (work evaluation) to run. Default : 10',
                        default=10)
    parser.add_argument('-maxh',
                        type=float, help='Terminate after this time. Default : 23.5 h',
                        default=23.5)
    parser.add_argument('--dry-run',
                        action='store_true', help='Print debug information')
    parser.add_argument('-MDRUN',
                        type=str, help='command for mdrun, we will use multidir, MPI is required. Default : mpirun -np 2 gmx_mpi mdrun',
                        default='mpirun -np 2 gmx_mpi mdrun')
    parser.add_argument('-GROMPP',
                        type=str, help='command for grompp, with additional flags. Default : gmx grompp',
                        default='gmx grompp')
    parser.add_argument('-tmp_folder',
                        type=str, help='Temporary folder. Point it to the local storage on the computing node to save IO.'
                                       ' You can also set this to /dev/shm if you have enough memory',)

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
                "dry_run"      : args.dry_run,
                "current_cycle": 0,
                "base_path"    : Path.cwd(),
                }
    if settings["dry_run"]:
        # set log to DEBUG
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
        for file_name in [settings["log"], settings["csv"]]:
            if file_name.exists():
                util.backup_if_exist(file_name)
        with open(settings["csv"], "w") as f:
            f.writelines(["Cycle,Work01(kJ/mol),Work10(kJ/mol),Acceptance_ratio,Accept\n"])
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
        w01, w10, p_accept= [float(i) for i in words[1:-1]]
        # check if cycle from csv matches the folder
        if cycle != settings["current_cycle"]:
            logging.info(f"Cycle from csv {cycle} does not match the folder {settings['current_cycle']}")
            exit(1)
        else:
            current_folder = Path(f"{settings['current_cycle']:06d}")
        if words[-1] == "A":
            swap = True
            s0 = current_folder / "1/ti.cpt", current_folder / "1/ti.gro"
            s1 = current_folder / "0/ti.cpt", current_folder / "0/ti.gro"
        elif words[-1] == "R":
            swap = False
            s0 = current_folder / "0/ti.cpt", current_folder / "0/ti.gro"
            s1 = current_folder / "1/ti.cpt", current_folder / "1/ti.gro"
        else:
            logging.info(f"Invalid csv line {l}, cannot proceed with restart.")
            exit(1)
        settings['current_cycle'] += 1
        current_folder = Path(f"{settings['current_cycle']:06d}")

    # prepare scratch folder, tmp_folder/0, tmp_folder/1
    if not settings["tmp_folder"].exists():
        settings["tmp_folder"].mkdir()
    for i in range(2):
        tmp_folder = settings["tmp_folder"] / str(i)
        if not tmp_folder.exists():
            tmp_folder.mkdir()
    # make sure the folder is empty
    for i in range(2):
        tmp_folder = settings["tmp_folder"] / str(i)
        for f in tmp_folder.iterdir():
            if f.is_file():
                f.unlink()
            else:
                shutil.rmtree(f)

    kBT = 8.314462618e-3 * settings["ref_t"] # kJ * K/mol
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
    logging.info(f"# Simulation settings End #########################################################################")




    # safety check, backup if needed

    # eq+NCMC cycle
    for i in range(settings["cycle"]):
        # check time
        if (time.time() - t0) / 3600 > settings["maxh"]:
            logging.info(f"Time limit reached, exit")
            break

        settings["current_folder"] = Path(f"{settings['current_cycle']:06d}")
        # eq
        logging.info(f"Cycle {settings['current_cycle']}, eq")
        # mkdir current_folder/0, current_folder/1
        for i in range(2):
            rep_dir = settings["current_folder"] / str(i)
            if not rep_dir.exists():
                rep_dir.mkdir(parents=True)

        if settings["current_cycle"]==0:
            run_eq_mdrun(settings)
        else:
            run_eq_grompp(settings, s0, s1)
            run_eq_mdrun(settings)

        # TI
        logging.info(f"Cycle {settings['current_cycle']}, TI")
        run_ti_grompp(settings)
        w01, w10 = run_ti_mdrun(settings)

        # swap attempt
        csv_line, s0, s1 =  swap_check(w01, w10, kBT, settings)
        with open(settings["csv"], "a") as f:
            f.write(csv_line)

        settings["current_cycle"] += 1
    logging.info(f"{settings['cycle']} cycle(s) have finished in {(time.time()-t0)/3600:.1f} h")
    logging.info(f"Cleaning up ...")
    shutil.rmtree(settings["tmp_folder"])
    logging.info(f"ALL Done")


