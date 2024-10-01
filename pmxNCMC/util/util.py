import os
import shutil
import logging

def backup_if_exist(f_name):
    """
    Do gromacs style backup. If md.xtc exist, backup it to #md.xtc.1#
    given f_name file should exist
    """
    for i in range(1,10000):
        bak = f"#{f_name}.{i}#"
        if not os.path.exists(bak):
            shutil.move(f_name, bak)
            logging.info(f"File {f_name} exists, backup to {bak}")
            return bak
    raise Exception(f"Cannot backup {f_name}")

def get_ref_T(mdp_file):
    """
    Get reference temperature from mdp file
    """
    with open(mdp_file) as f:
        for line in f:
            if "ref_t" in line or "ref-t" in line:
                l = line.split(";")[0]
                words = l.split()
                if len(words) != 3:
                    logging.info(f"There are more than 1 ref_t in {line.rstrip()}, return the first one")
                return float(words[2])
    raise Exception(f"Cannot find ref_t in {mdp_file}")


def mdp_check_TI(mdp_file):
    """
    Check if the mdp file has the correct init-lambda, delta-lambda, nsteps,
    if init-lambda = 1. delta-lambda * nsteps = -1
    if init-lambda = 0. delta-lambda * nsteps = 1
    :param mdp_lines: 
    :return: T/F
    """
    with open(mdp_file) as f:
        mdp_lines = f.readlines()
    init_lambda = None
    delta_lambda = None
    nsteps = None
    for line in mdp_lines:
        l = line.split(";")[0]
        if "init-lambda" in l or "init_lambda" in l:
            init_lambda = float(l.split()[2])
        if "delta-lambda" in l or "delta_lambda" in l:
            delta_lambda = float(l.split()[2])
        if "nsteps" in l:
            nsteps = int(l.split()[2])
    if init_lambda is None or delta_lambda is None or nsteps is None:
        logging.info(f"{init_lambda} {delta_lambda} {nsteps}")
        logging.info(f"Cannot find init-lambda, delta-lambda, or nsteps in {mdp_file}")
        return False
    if init_lambda == 1 and delta_lambda * nsteps == -1:
        return True
    if init_lambda == 0 and delta_lambda * nsteps == 1:
        return True
    logging.info(f"MDP file {mdp_file} does not have the correct init-lambda, delta-lambda, nsteps")
    logging.info(f"init-lambda={init_lambda}, delta-lambda * nsteps = {delta_lambda * nsteps}")
    return False