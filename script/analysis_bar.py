#!/usr/bin/env python3

import pmxNCMC
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymbar

def gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y

def data2gauss(data):
    '''Takes a one dimensional array and fits a Gaussian.

    Returns
    -------
    float
        mean of the distribution.
    float
        standard deviation of the distribution.
    float
        height of the curve's peak.
    '''
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A

def plot_work_dist(wf, wr, fname="Wdist.png", nbins=20, dG=None, dGerr=None, units=None, kBT=None, kBT_gmx=None):
    """
    Plot the work distribution
    :param wf:
    :param wr:
    :param fname:
    :param nbins:
    :param dG:    kB*T unit
    :param dGerr: kB*T unit
    :param units:
    :param kBT:
    :return:
    """
    w_f = wf.values/kBT_gmx*kBT  # to target unit
    w_r = wr.values/kBT_gmx*kBT
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))
    x1 = np.arange(len(w_f))
    sm1 = np.convolve(w_f, np.ones(11) / 11, mode='valid')
    x2 = np.arange(len(w_r))
    sm2 = np.convolve(w_r, np.ones(11) / 11, mode='valid')

    ax1.plot(x1,       w_f, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
    ax1.plot(x1[5:-5], sm1, 'g-', linewidth=3)
    ax1.plot(x2,       w_r, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
    ax1.plot(x2[5:-5], sm2, 'b-', linewidth=3)
    ax1.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    ax1.set_ylabel(f'W {units}', fontsize=20)
    ax1.set_xlabel(r'# Snapshot', fontsize=20)
    ax1.grid(lw=2)
    ax1.set_xlim(0, max(len(x1), len(x2)) + 1)

    bins = np.linspace(min(w_f.min(), w_r.min()), max(w_f.max(), w_r.max()), nbins)
    ax2.hist(w_f, bins=bins, orientation='horizontal', facecolor='green',
             alpha=.75, density=True)
    ax2.hist(w_r, bins=bins, orientation='horizontal', facecolor='blue',
             alpha=.75, density=True)

    x = np.linspace(min(w_f.min(), w_r.min()), max(w_f.max(), w_r.max()), 1000)
    mf, devf, Af = data2gauss(w_f)
    mb, devb, Ab = data2gauss(w_r)
    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)
    ax2.plot(y1, x, 'g--', linewidth=2)
    ax2.plot(y2, x, 'b--', linewidth=2)
    size = max([max(y1), max(y2)])
    res_x = [dG*kBT, dG*kBT]
    res_y = [0, size*1.2]
    if dG is not None and dGerr is not None:
        ax2.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f $\pm$ %.2f %s' % (dG*kBT, dGerr*kBT, units))
        ax2.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    fig.savefig(fname, dpi=300)

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
                        choices=["kJ", "kj", "kJ/mol", "kj/mol",
                                 "kcal", "kCal", "kcal/mol", "kCal/mol"])
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
        plot_work_dist(work01, -1 * work10, fname=args.w, nbins=20, dG=dG, dGerr=dGe, units=unit, kBT=kBT, kBT_gmx=kBT_gmx)

