import os
import shutil
import logging
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
logging.getLogger("pymbar").setLevel(logging.ERROR) # Suppress logging in pymbar
import pymbar

kB_kj_mol   = 8.314462618e-3
kB_kcal_mol = 1.987204259e-3

def integrate_work(ti_xvg):
    """
    Integrate the work from a ti.xvg file. Only works for linear lambda 0->1 or 1->0
    :param ti_xvg: str, path to ti.xvg file
    :return: integrated work from 0 to 1. kJ/mol, the same as gromacs
    """
    dh_dl = np.loadtxt(ti_xvg, comments=["@", "#"])
    lam = np.linspace(0, 1, len(dh_dl))
    return simpson(dh_dl[:, 1], x=lam)


def backup_if_exist(f_name):
    """
    Do gromacs style backup. If md.xtc exist, backup it to #md.xtc.1#
    given f_name file should exist
    """
    for i in range(1,10000):
        bak = f"#{f_name}.{i}#"
        if not os.path.exists(bak):
            shutil.move(f_name, bak)
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
                t_words = l.split("=")[1].split()
                if len(t_words) != 1:
                    logging.info(f"There are more than 1 ref_t in {line.rstrip()}, return the first one")
                return float(t_words[0])
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
            init_lambda = float(l.split("=")[1])
        if "delta-lambda" in l or "delta_lambda" in l:
            delta_lambda = float(l.split("=")[1])
        if "nsteps" in l:
            nsteps = int(l.split("=")[1])
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

def free_E_bar(work01, work10):
    """
    Given the work in the two direction, estimate the free energy using BAR
    everything should be in kBT unit
    Will try pymbar 3 and 4.
    :param work01:
    :param work10:
    :return: dG, dG_error
    """
    try:
        # pymbar 3
        dG, dGe = pymbar.BAR(work01, work10)
    except AttributeError as e:
        # pymbar 4
        res = pymbar.other_estimators.bar(work01, work10)
        dG = res["Delta_f"]
        dGe = res["dDelta_f"]
    return dG, dGe

# for plotting
def gauss_func(A, mean, dev, x):
    """Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function"""
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y

def data2gauss(data):
    """
    Given data, return the mean, deviation and amplitude of the Gaussian
    :param data: 1D array
    :return: mean, deviation, amplitude
    """
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A

def plot_work_dist_line(ax, w_f, w_r,  smooth_window=11):
    """
    Plot the work distribution in a line plot
    :param ax: matplotlib axis
    :param w_f: Forward work in the target unit
    :param w_r: Reverse work in the target unit
    :param smooth_window: size of the window in moving window average (smoothing).
        If there is less data, smoothed line will not be plotted
    :return: None
    """
    x1 = np.arange(len(w_f))
    x2 = np.arange(len(w_r))
    if len(x1) > smooth_window+1:
        ax.plot(x1, w_f, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
        sm1 = np.convolve(w_f, np.ones(11) / 11, mode='valid')
        ax.plot(x1[5:-5], sm1, 'g-', linewidth=3)
    else:
        ax.plot(x1, w_f, 'g-', linewidth=3, label="Forward (0->1)")
    if len(x2) > smooth_window+1:
        ax.plot(x2, w_r, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
        sm2 = np.convolve(w_r, np.ones(11) / 11, mode='valid')
        ax.plot(x2[5:-5], sm2, 'b-', linewidth=3)
    else:
        ax.plot(x2, w_r, 'b-', linewidth=3, label="Backward (1->0)")
    ax.set_xlim(0, max(len(x1), len(x2)) + 1)

def plot_work_dist_hist(ax, w_f, w_r, bins=20):
    """
    Plot the work distribution in a histogram in the horizontal orientation, so that it can be compared side by side with the line plot
    :param ax: matplotlib axis
    :param w_f: Forward work in the target unit
    :param w_r: Reverse work in the target unit
    :param bins: for the histogram
    :return:
    """
    ax.hist(w_f, bins=bins, orientation='horizontal', facecolor='green',
             alpha=.75, density=True)
    ax.hist(w_r, bins=bins, orientation='horizontal', facecolor='blue',
             alpha=.75, density=True)
    x = np.linspace(min(w_f.min(), w_r.min()), max(w_f.max(), w_r.max()), 1000)
    mf, devf, Af = data2gauss(w_f)
    mb, devb, Ab = data2gauss(w_r)
    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)
    ax.plot(y1, x, 'g--', linewidth=2)
    ax.plot(y2, x, 'b--', linewidth=2)

def plot_work_dist(wf, wr, kBT_in=None, fname="Wdist.png", nbins=20, dG=None, dGerr=None, units=None, kBT=None, ):
    """
    Plot the work distribution
    In ax1, line plot of forward and backward work
    In ax2, histogram of two work distribution in the horizontal orientation
    :param wf: Forward work. This should use the same unit as kBT_in
    :param wr: Backward work. This should use the same unit as kBT_in
    :param kBT_in: kB*T in the unit of input work
    :param fname: file name to save the plot
    :param nbins: number of bins for the histogram
    :param dG:    kB*T unit
    :param dGerr: kB*T unit
    :param units: (str) name of tje unit in the output
    :param kBT:  kB*T in the unit of output
    For unit conversion, we assume wf/kBT_in * kBT gives the target unit

    :return:
    """
    w_f = wf.values / kBT_in * kBT  # to target unit
    w_r = wr.values / kBT_in * kBT
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))
    plot_work_dist_line(ax1, w_f, w_r)
    ax1.legend(loc='upper center', prop={'size': 12})
    ax1.set_ylabel(f'W [{units}]', fontsize=20)
    ax1.set_xlabel(r'# Snapshot', fontsize=20)
    ax1.grid(lw=2)


    bins = np.linspace(min(w_f.min(), w_r.min()), max(w_f.max(), w_r.max()), nbins)
    plot_work_dist_hist(ax2, w_f, w_r, bins)
    res_x = [dG*kBT, dG*kBT]
    res_y = [0, ax2.get_xlim()[1]*0.95]
    if dG is not None and dGerr is not None:
        ax2.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f $\pm$ %.2f %s' % (dG*kBT, dGerr*kBT, units))
        ax2.legend(loc='upper center', prop={'size': 12})
    ax2.set_ylim(ax1.get_ylim())
    fig.savefig(fname, dpi=300)