
#
#   Plot and export filter or smoother position and velocity covariance
#   (JSON input version)
#
#   Reads a GMAT filter or smoother JSON output file.  CovarianceVNB is taken
#   directly from the file when the run used DataFileStyle = 'Verbose'; otherwise
#   it is computed from the Cartesian covariance and position/velocity state.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the JSON file.
#
#   usage: plot_covariance_json.py [-h] [--no_plot] [--no_csv] jsonfile
#
#   jsonfile    (required) full path to filter or smoother JSON output file
#   --no_plot   (optional) don't show plot or create plot PDF files
#   --no_csv    (optional) don't create the CSV data files
#   -h, --help  (optional) show help message
#

import json
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, sys, os

from xyz_to_vnb import xyz_to_vnb


def parse_epoch(s):
    """Parse a GMAT UTC epoch string like '10 Jun 2010 00:00:00.000'."""
    return datetime.strptime(s, "%d %b %Y %H:%M:%S.%f")


def load(jsonfile):

    with open(jsonfile) as f:
        data = json.load(f)

    if 'FilterUpdates' in data:
        runtype = 'Filter'
        updates = data['FilterUpdates']
    elif 'SmootherUpdates' in data:
        runtype = 'Smoother'
        updates = data['SmootherUpdates']
    else:
        print('Error: Cannot find FilterUpdates or SmootherUpdates in', jsonfile)
        sys.exit()

    epochs = []
    cov_pos_v, cov_pos_n, cov_pos_b = [], [], []
    cov_vel_v, cov_vel_n, cov_vel_b = [], [], []

    for u in updates:

        cov_cart = np.array(u['Covariance'])

        if 'CovarianceVNB' in u:
            cov_vnb = np.array(u['CovarianceVNB'])
        else:
            state = np.array(u['State'])
            M = xyz_to_vnb(state)
            cov_vnb = M @ cov_cart[0:6, 0:6] @ M.T

        epochs.append(parse_epoch(u['Epoch']))
        cov_pos_v.append(cov_vnb[0, 0])
        cov_pos_n.append(cov_vnb[1, 1])
        cov_pos_b.append(cov_vnb[2, 2])
        cov_vel_v.append(cov_vnb[3, 3])
        cov_vel_n.append(cov_vnb[4, 4])
        cov_vel_b.append(cov_vnb[5, 5])

    return epochs, runtype, \
        [cov_pos_v, cov_pos_n, cov_pos_b, cov_vel_v, cov_vel_n, cov_vel_b]


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    t, runtype, \
        [cov_pos_v, cov_pos_n, cov_pos_b, cov_vel_v, cov_vel_n, cov_vel_b] = data

    #
    #   Plot position covariance
    #

    f_pos = plt.figure()
    f_pos.suptitle(runtype + ' 1-Sigma Position Uncertainty')

    plt.plot_date(t, np.sqrt(cov_pos_v) * 1e3, label='Sigma-V',
        xdate=True, linestyle='solid', fmt='C0')
    plt.plot_date(t, np.sqrt(cov_pos_n) * 1e3, label='Sigma-N',
        xdate=True, linestyle='solid', fmt='C1')
    plt.plot_date(t, np.sqrt(cov_pos_b) * 1e3, label='Sigma-B',
        xdate=True, linestyle='solid', fmt='C2')

    plt.legend(loc='lower center', ncol=3)
    plt.xlabel('Time UTC')
    plt.ylabel('meters')

    f_pos.autofmt_xdate()

    #
    #   Plot velocity covariance
    #

    f_vel = plt.figure()
    f_vel.suptitle(runtype + ' 1-Sigma Velocity Uncertainty')

    plt.plot_date(t, np.sqrt(cov_vel_v) * 1e5, label='Sigma-V',
        xdate=True, linestyle='solid', fmt='C0')
    plt.plot_date(t, np.sqrt(cov_vel_n) * 1e5, label='Sigma-N',
        xdate=True, linestyle='solid', fmt='C1')
    plt.plot_date(t, np.sqrt(cov_vel_b) * 1e5, label='Sigma-B',
        xdate=True, linestyle='solid', fmt='C2')

    plt.legend(loc='lower center', ncol=3)
    plt.xlabel('Time UTC')
    plt.ylabel('cm/sec')

    f_vel.autofmt_xdate()

    if save_plot:

        outfile = os.path.join(outdir, runtype + '_position_covariance.pdf')
        f_pos.savefig(outfile, bbox_inches='tight')

        outfile = os.path.join(outdir, runtype + '_velocity_covariance.pdf')
        f_vel.savefig(outfile, bbox_inches='tight')

    if save_csv:

        outfile = os.path.join(outdir, runtype + '_covariance.csv')
        with open(outfile, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Time',
                'Sigma-Pos-V (km)', 'Sigma-Pos-N (km)', 'Sigma-Pos-B (km)',
                'Sigma-Vel-V (km/sec)', 'Sigma-Vel-N (km/sec)', 'Sigma-Vel-B (km/sec)'])

            outarray = np.array([t,
                np.sqrt(cov_pos_v), np.sqrt(cov_pos_n), np.sqrt(cov_pos_b),
                np.sqrt(cov_vel_v), np.sqrt(cov_vel_n), np.sqrt(cov_vel_b)], dtype=object)

            for row in outarray.transpose():
                csv_writer.writerow(row)

    if show_plot:
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('jsonfile',
        help='Name of GMAT filter or smoother JSON output file')
    parser.add_argument('--no_plot',
        help="Don't show or save plot", action='store_false')
    parser.add_argument('--no_csv',
        help="Don't save data to CSV file", action='store_false')

    args = parser.parse_args()
    outdir, _ = os.path.split(args.jsonfile)

    data = load(args.jsonfile)

    render(data, outdir, show_plot=args.no_plot,
        save_plot=args.no_plot, save_csv=args.no_csv)
