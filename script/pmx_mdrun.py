#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path
import time
import logging

import pmxNCMC

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     f"""Version {pmxNCMC.__version__}. This is the IO based implementation for doing 
                                     expanded ensemble NCMC with only lambda 0 and 1""")
    parser.add_argument('-p',
                        type=str, help='Topology file',
                        default='topol.top')
    parser.add_argument('-log',
                        type=str, help='Log file',
                        default='md.log')
    parser.add_argument('-csv',
                        type=str, help='CSV file for saving work values',
                        default='md.csv')
    parser.add_argument('-mdp_folder',
                        type=str,
                        help='Folder that contains eq0.mdp, eq1.mdp, ti0.mdp, ti1.mdp. All 4 files are required')
    parser.add_argument('-folder_start',
                        type=str, help='Folder to start the simulation.',
                        default='00000')
    parser.add_argument('-cycle',
                        type=int, help='Number of cycles (work evaluation) to run',
                        default=10)
    parser.add_argument('-maxh',
                        type=float, help='Terminate after 0.99 times this time')
    parser.add_argument('--dry-run',
                        action='store_true', help='Test 1 cycle')
    parser.add_argument('-MDRUN',
                        type=str, help='command for mdrun, we will use multidir, MPI is required',
                        default='mpirun -np 2 gmx_mpi mdrun')
    parser.add_argument('-GROMPP',
                        type=str, help='command for grompp, with additional flags',
                        default='gmx_mpi grompp')
    parser.add_argument('-tmp_folder',
                        type=str, help='Temporary folder. Point it to the local storage on the computing node to save IO.'
                                       ' You can also set this to /dev/shm if you have enough memory',
                        default="./")

    args = parser.parse_args()
    settings = {"top"          : args.p,
                "log"          : args.log,
                "csv"          : args.csv,
                "mdp_folder"   : Path(args.mdp_folder),
                "folder_start" : Path(args.folder_start),
                "cycle"        : args.cycle,
                "maxh"         : args.maxh,
                "MDRUN"        : args.MDRUN,
                "tmp_folder"   : args.tmp_folder,
                "dry_run"      : args.dry_run}

    # safety check, backup if needed

    # eq+NCMC cycle
    for i in range(settings["cycle"]):
        # eq
        # NCMC
        pass
