
#
#   Plot and export filter/smoother consistency tests
#   (JSON input version)
#
#   Generates filter/smoother consistency plots for position, velocity, and all
#   solve-fors in the run.  Filter and smoother updates are matched by epoch.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the filter JSON file.
#
#   usage: filter_smoother_consistency_json.py [-h] [--no_plot] [--no_csv]
#                                     filter_jsonfile smoother_jsonfile
#
#   filter_jsonfile    (required) full path to filter JSON output file
#   smoother_jsonfile  (required) full path to smoother JSON output file
#   --no_plot          (optional) don't show plot or create plot PDF files
#   --no_csv           (optional) don't create the CSV data files
#   -h, --help         (optional) show help message
#

import json
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, os, sys

from xyz_to_vnb import xyz_to_vnb


def parse_epoch(s):
    """Parse a GMAT UTC epoch string like '10 Jun 2010 00:00:00.000'."""
    return datetime.strptime(s, "%d %b %Y %H:%M:%S.%f")


def load(filter_jsonfile, smoother_jsonfile):

    with open(filter_jsonfile) as f:
        fildata = json.load(f)
    with open(smoother_jsonfile) as f:
        smodata = json.load(f)

    if 'FilterUpdates' not in fildata:
        print('Error: FilterUpdates not found in', filter_jsonfile)
        sys.exit()
    if 'SmootherUpdates' not in smodata:
        print('Error: SmootherUpdates not found in', smoother_jsonfile)
        sys.exit()

    state_names = fildata['CartesianStateNames']
    nstates = len(state_names)

    # Build smoother lookup keyed by epoch string
    smo_by_epoch = {u['Epoch']: u for u in smodata['SmootherUpdates']}

    epochs   = []
    fil_states = []
    fil_covs   = []
    smo_states = []
    smo_covs   = []

    for u in fildata['FilterUpdates']:
        ep = u['Epoch']
        if ep not in smo_by_epoch:
            continue
        su = smo_by_epoch[ep]
        epochs.append(parse_epoch(ep))
        fil_states.append(np.array(u['State']))
        fil_covs.append(np.array(u['Covariance']))
        smo_states.append(np.array(su['State']))
        smo_covs.append(np.array(su['Covariance']))

    if not epochs:
        print('Error: No matching epochs found between filter and smoother.')
        sys.exit()

    # Omit last epoch (square-root of near-zero difference is unstable there)
    n = len(epochs) - 1
    epochs     = epochs[:n]
    fil_states = fil_states[:n]
    fil_covs   = fil_covs[:n]
    smo_states = smo_states[:n]
    smo_covs   = smo_covs[:n]

    #
    #   Position/velocity consistency in VNB frame
    #

    posvel_consistency = np.zeros((6, n))

    for i in range(n):
        rv_fil = fil_states[i][0:6]
        rv_smo = smo_states[i][0:6]
        M = xyz_to_vnb(rv_fil)

        rv_fil_vnb = M @ rv_fil
        rv_smo_vnb = M @ rv_smo

        fil_cov_vnb = M @ fil_covs[i][0:6, 0:6] @ M.T
        smo_cov_vnb = M @ smo_covs[i][0:6, 0:6] @ M.T

        delta   = rv_fil_vnb - rv_smo_vnb
        sigma   = np.sqrt(np.diag(fil_cov_vnb) - np.diag(smo_cov_vnb))
        posvel_consistency[:, i] = np.divide(delta, sigma)

    #
    #   Consistency for any additional solve-fors
    #

    param_names = state_names[6:]
    param_consistency = np.zeros((len(param_names), n))

    for k, idx in enumerate(range(6, nstates)):
        for i in range(n):
            delta = fil_states[i][idx] - smo_states[i][idx]
            sigma = np.sqrt(fil_covs[i][idx, idx] - smo_covs[i][idx, idx])
            param_consistency[k, i] = delta / sigma

    return epochs, posvel_consistency, param_names, param_consistency


def _sigma_lines(t, sigma=3):
    return [+sigma] * len(t), [-sigma] * len(t)


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    t, posvel_consistency, param_names, param_consistency = data

    plus3, minus3 = _sigma_lines(t)

    #
    #   Position consistency
    #

    pos_c = posvel_consistency[0:3, :]

    f = plt.figure()
    f.suptitle('Filter-Smoother Position Consistency')

    plt.plot_date(t, pos_c[0, :], label='Consistency-V',
        xdate=True, linestyle='solid', fmt='C0')
    plt.plot_date(t, pos_c[1, :], label='Consistency-N',
        xdate=True, linestyle='solid', fmt='C1')
    plt.plot_date(t, pos_c[2, :], label='Consistency-B',
        xdate=True, linestyle='solid', fmt='C2')

    plt.plot_date(t,  plus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')
    plt.plot_date(t, minus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')

    plt.legend(loc='lower center', ncol=4, mode='expand')
    plt.xlabel('Time UTC')
    plt.ylabel('Unitless')
    f.autofmt_xdate()

    if save_plot:
        f.savefig(os.path.join(outdir, 'fs_position_consistency.pdf'), bbox_inches='tight')

    if save_csv:
        outfile = os.path.join(outdir, 'fs_position_consistency.csv')
        with open(outfile, 'w', newline='') as csvfile:
            w = csv.writer(csvfile)
            w.writerow(['Time', 'Consistency-V', 'Consistency-N', 'Consistency-B'])
            for row in np.array([t, pos_c[0], pos_c[1], pos_c[2]], dtype=object).T:
                w.writerow(row)

    #
    #   Velocity consistency
    #

    vel_c = posvel_consistency[3:6, :]

    f = plt.figure()
    f.suptitle('Filter-Smoother Velocity Consistency')

    plt.plot_date(t, vel_c[0, :], label='Consistency-V',
        xdate=True, linestyle='solid', fmt='C0')
    plt.plot_date(t, vel_c[1, :], label='Consistency-N',
        xdate=True, linestyle='solid', fmt='C1')
    plt.plot_date(t, vel_c[2, :], label='Consistency-B',
        xdate=True, linestyle='solid', fmt='C2')

    plt.plot_date(t,  plus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')
    plt.plot_date(t, minus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')

    plt.legend(loc='lower center', ncol=4, mode='expand')
    plt.xlabel('Time UTC')
    plt.ylabel('Unitless')
    f.autofmt_xdate()

    if save_plot:
        f.savefig(os.path.join(outdir, 'fs_velocity_consistency.pdf'), bbox_inches='tight')

    if save_csv:
        outfile = os.path.join(outdir, 'fs_velocity_consistency.csv')
        with open(outfile, 'w', newline='') as csvfile:
            w = csv.writer(csvfile)
            w.writerow(['Time', 'Consistency-V', 'Consistency-N', 'Consistency-B'])
            for row in np.array([t, vel_c[0], vel_c[1], vel_c[2]], dtype=object).T:
                w.writerow(row)

    #
    #   Any additional solve-fors
    #

    for i, param in enumerate(param_names):

        f = plt.figure()
        f.suptitle('Filter-Smoother ' + param + ' Consistency')

        plt.plot_date(t, param_consistency[i, :], label=param,
            xdate=True, linestyle='solid', fmt='C0')

        plt.plot_date(t,  plus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')
        plt.plot_date(t, minus3, linestyle=(0, (5, 5)), alpha=0.5, linewidth=0.7, fmt='k')

        plt.legend(loc='lower center', ncol=4, mode='expand')
        plt.xlabel('Time UTC')
        plt.ylabel('Unitless')
        f.autofmt_xdate()

        if save_plot:
            f.savefig(os.path.join(outdir, 'fs_' + param + '_consistency.pdf'),
                bbox_inches='tight')

        if save_csv:
            outfile = os.path.join(outdir, 'fs_' + param + '_consistency.csv')
            with open(outfile, 'w', newline='') as csvfile:
                w = csv.writer(csvfile)
                w.writerow(['Time', param + ' Consistency'])
                for row in np.array([t, param_consistency[i]], dtype=object).T:
                    w.writerow(row)

    if show_plot:
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('filter_jsonfile',
        help='Name of GMAT filter JSON output file')
    parser.add_argument('smoother_jsonfile',
        help='Name of GMAT smoother JSON output file')
    parser.add_argument('--no_plot',
        help="Don't show or save plot", action='store_false')
    parser.add_argument('--no_csv',
        help="Don't save data to CSV file", action='store_false')

    args = parser.parse_args()
    outdir, _ = os.path.split(args.filter_jsonfile)

    data = load(args.filter_jsonfile, args.smoother_jsonfile)

    render(data, outdir, show_plot=args.no_plot,
        save_plot=args.no_plot, save_csv=args.no_csv)
