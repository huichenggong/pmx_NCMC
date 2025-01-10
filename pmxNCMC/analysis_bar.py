#!/usr/bin/env python3
import logging
import sys
import argparse

import numpy as np
import pandas as pd

import pmxNCMC

logging.getLogger("pymbar").setLevel(logging.ERROR) # Suppress logging in pymbar
import pymbar


def main():
    parser = argparse.ArgumentParser(
        description=f"""Version {pmxNCMC.__version__}. Read the csv file from pmx_mdrun.py and estimate the free energy difference using BAR.
                            Sample usage: analysis_bar.py -csv md.csv -w wplot.png --unit kcal 
                            # Work and temperature will be read from md.csv. The work for each snapshot will be plotted in wplot.png. The unit is kcal/mol.""")
    parser.add_argument("-csv",
                        metavar='csv input',
                        type=str, help="csv file from pmx_mdrun.py. The 2nd, 3rd columns are the work in kJ.mol. "
                                       "Default : md.csv",
                        default="md.csv")
    parser.add_argument("-t",
                        metavar='temperature',
                        type=float, help='Temperature in K, If not given, try to read from the heading of the csv')
    parser.add_argument("-oA",
                        metavar='integA.dat',
                        type=str, help='work output A->B. File to save all the work value (for pmx). The unis is kJ/mol. ')
    parser.add_argument("-oB",
                        metavar='integB.dat',
                        type=str, help='work output B->A. File to save all the work value (for pmx). The unis is kJ/mol. ')
    parser.add_argument("--unit",
                        metavar='unit',
                        type=lambda x: x.lower(),
                        help='Unit for the print out and plotting. Default is "kJ/mol"',
                        default="kj",
                        choices=["kj", "kj/mol", "kcal", "kcal/mol"])
    parser.add_argument("-w",
                        metavar='plot',
                        type=str, help='Work distribution plot. Default is "Wplot.png"',
                        default="Wplot.png")
    parser.add_argument("-b",
                        metavar="begin",
                        type=int, help="The beginning of the data. Default is 0",
                        default=0)
    parser.add_argument("-e",
                        metavar="end",
                        type=int, help="The end of the data. Default is the end of the data")

    args = parser.parse_args()

    print(f"pmx_mdrun version {pmxNCMC.__version__}")
    print(f"Command line arguments:")
    print(f"  {' '.join(sys.argv)}")
    print(f"Reading {args.csv}")
    df = pd.read_csv(args.csv)

    # get temperature
    if len(df.columns[-1]) > 7:
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
    kBT_gmx = pmxNCMC.util.kB_kj_mol * temperature  # kJ/mol

    # unit conversion
    if args.unit.lower() in ["kj", "kj/mol"]:
        unit = "kJ/mol"
        kB = pmxNCMC.util.kB_kj_mol
    elif args.unit.lower() in ["kcal", "kcal/mol"]:
        unit = "kcal/mol"
        kB = pmxNCMC.util.kB_kcal_mol
    else:
        raise ValueError(f"Units should be kJ/mol or kcal/mol. {args.unit} is not supported")
    kBT = kB * temperature
    print(f"The selected unit is {unit}.")
    print(f"kB  = {kB:.5e} {unit}/T")
    print(f"kBT = {kBT:.2f} {unit}")

    # BAR estimation
    work01 = df[df.columns[1]]  # 0->1, kJ/mol
    work10 = df[df.columns[2]]
    if args.e:
        work01 = work01[:args.e]
        work10 = work10[:args.e]
    if args.b:
        work01 = work01[args.b:]
        work10 = work10[args.b:]

    print( "Number of work in 0->1 , 1->0         : ", len(work01), len(work10))
    print(f"Number of Accepted/Attempted exchange : { sum( df[df.columns[-1]] == 'A' ) } / {len(df)}")
    dG, dGe = pmxNCMC.util.free_E_bar(work01/kBT_gmx, work10/kBT_gmx)
    print(f"DeltaG = {dG * kBT:.2f} +- {dGe * kBT:.2f} {unit}")
    try:
        overlap = pymbar.bar_overlap(work01 / kBT, work10 / kBT)
        print(f"Overlap BAR: {overlap:.2f}")
    except AssertionError:
        print(f"Overlap BAR: 0 (pymbar failed to calculate the overlap)")

    # save work for pmx
    if args.oA:
        work01.to_csv(args.oA, sep=' ', header=False)
    if args.oB:
        (work10 * -1).to_csv(args.oB, sep=' ', header=False)

    # plot
    if args.w not in ["None", "none"]:
        pmxNCMC.util.plot_work_dist(work01, -1 * work10, kBT_in=kBT_gmx,
                                    fname=args.w, nbins=20, dG=dG, dGerr=dGe, units=unit, kBT=kBT)


if __name__ == "__main__":
    main()