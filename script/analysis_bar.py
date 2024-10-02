#!/usr/bin/env python3

import pmxNCMC
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
                        type=str, help="csv file from pmx_mdrun.py Default : md.csv", required=True,
                        default="md.csv")
    parser.add_argument("-t",
                        metavar='temperature',
                        type=float, help='Temperature in K, Default is 298.15 K', default=298.15)
    parser.add_argument("-oA",
                        metavar='work output 0->1',
                        type=str, help='File.dat where to save all the work value (for pmx). Default is "integA.dat"')
    parser.add_argument("-oB",
                        metavar='work output 1->0',
                        type=str, help='File.dat where to save all the work value (for pmx). Default is "integB.dat"')
    parser.add_argument("--unit",
                        metavar='unit',
                        type=str, help='Unit for the work. Default is "kJ/mol"',
                        default="kJ",
                        choices=["kJ", "kj", "KJ", "Kj",
                                 "kJ/mol", "kj/mol", "KJ/mol", "Kj/mol",
                                 "kcal", "kCal", "Kcal", "KCal",
                                 "kcal/mol", "kCal/mol", "Kcal/mol", "KCal/mol"])
    parser.add_argument("-w",
                        metavar='plot',
                        type=str, help='Work distribution plot. Default is "wplot.png"')

    args = parser.parse_args()

    if args.unit.lower() in ["kj", "kj/mol"]:
        unit = "kJ/mol"
        kBT = 8.314462618e-3 * args.t
    elif args.unit.lower() in ["kcal", "kcal/mol"]:
        unit = "kcal/mol"
        kBT = 1.987204259e-3 * args.t
    else:
        raise ValueError(f"Units should be kJ/mol or kcal/mol. {args.unit} is not supported")


    kBT_gmx = 8.314462618e-3 * args.t  # kJ/mol
    df = pd.read_csv(args.csv)
    heading = list(df.keys())
    work01 = df[heading[1]]
    work10 = df[heading[2]]
    try:
        # pymbar 3
        dG, dGe = pymbar.BAR(work01 / kBT_gmx, work10 / kBT_gmx)
    except AttributeError as e:
        # pymbar 4
        res = pymbar.other_estimators.bar(work01 / kBT_gmx, work10 / kBT_gmx)
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

