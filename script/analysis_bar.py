#!/usr/bin/env python3

import pmxNCMC
import sys
import argparse
import numpy as np
import pandas as pd
import logging
logging.getLogger("pymbar").setLevel(logging.ERROR) # Suppress logging in pymbar
import pymbar



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""Version {pmxNCMC.__version__}.""")
    parser.add_argument("-csv",
                        metavar='csv input',
                        type=str, help="csv file from pmx_mdrun.py. The 2nd, 3rd columns are the work in kJ.mol. "
                                       "Default : md.csv",
                        default="md.csv")
    parser.add_argument("-t",
                        metavar='temperature',
                        type=float, help='Temperature in K, If not given, try to read from the heading of the csv')
    parser.add_argument("-oA",
                        metavar='work output 0->1',
                        type=str, help='File.dat where to save all the work value (for pmx). The unit is kJ/mol. '
                                       'Default is "integA.dat"')
    parser.add_argument("-oB",
                        metavar='work output 1->0',
                        type=str, help='File.dat where to save all the work value (for pmx). The unis is kJ/mol. '
                                       'Default is "integB.dat"')
    parser.add_argument("--unit",
                        metavar='unit',
                        type=str, help='Unit for the print out and plotting. Default is "kJ/mol"',
                        default="kJ",
                        choices=["kJ", "kj", "KJ", "Kj",
                                 "kJ/mol", "kj/mol", "KJ/mol", "Kj/mol",
                                 "kcal", "kCal", "Kcal", "KCal",
                                 "kcal/mol", "kCal/mol", "Kcal/mol", "KCal/mol"])
    parser.add_argument("-w",
                        metavar='plot',
                        type=str, help='Work distribution plot. Default is "wplot.png"')

    args = parser.parse_args()




    print(f"Command line arguments: {' '.join(sys.argv)}")
    print(f"Reading {args.csv}")
    df = pd.read_csv(args.csv)
    if len(df.columns[-1])>7:
        words = df.columns[-1].split("_")
        if len(words) == 2:
            temperature = float(df.columns[-1].split("_")[-1])
        else:
            raise ValueError(f"Cannot read temperature from the last column name: {df.columns[-1]}")
        if args.t:
            np.allclose(temperature, args.t, atol=1e-4)
            print(f"Temperature is the same in the -t and {args.csv}: {temperature} K")
        else:
            print(f"Temperature is read from the csv file: {temperature} K")
    else:
        if not args.t:
            raise ValueError(f"Temperature is not found in the csv file, please give it by -t")
        else:
            temperature = args.t
            print(f"Temperature is given by argument -t : {args.t} K")
    kBT_gmx = 8.314462618e-3 * temperature  # kJ/mol
    # unit conversion
    if args.unit.lower() in ["kj", "kj/mol"]:
        unit = "kJ/mol"
        kB = 8.314462618e-3
    elif args.unit.lower() in ["kcal", "kcal/mol"]:
        unit = "kcal/mol"
        kB = 1.987204259e-3
    else:
        raise ValueError(f"Units should be kJ/mol or kcal/mol. {args.unit} is not supported")
    kBT = kB * temperature
    print(f"The selected unit is {unit}.")
    print(f"kB  = {kB:.5e} {unit}/T")
    print(f"kBT = {kBT:.2f} {unit}")

    work01 = df[df.columns[1]] # 0->1, kJ/mol
    work10 = df[df.columns[2]]
    try:
        # pymbar 3
        dG, dGe = pymbar.BAR(work01 / kBT_gmx, work10 / kBT_gmx)
        print("pymbar 3 is used. >>> pymbar.BAR(w_f, w_r)")
    except AttributeError as e:
        # pymbar 4
        res = pymbar.other_estimators.bar(work01 / kBT_gmx, work10 / kBT_gmx)
        print("pymbar 4 is used. >>> pymbar.other_estimators.bar(w_f, w_r)")
        dG = res["Delta_f"]
        dGe = res["dDelta_f"]

    print(f"DeltaG = {dG*kBT:.2f} +- {dGe*kBT:.2f} {unit}")

    if args.oA:
        work01.to_csv(args.oA, sep=' ', header=False)
    if args.oB:
        (work10*-1).to_csv(args.oB, sep=' ', header=False)

    if args.w:
        pmxNCMC.util.plot_work_dist(work01, -1 * work10, kBT_in=kBT_gmx,
                                    fname=args.w, nbins=20, dG=dG, dGerr=dGe, units=unit, kBT=kBT)

