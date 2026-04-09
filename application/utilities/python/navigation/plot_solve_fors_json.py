
#
#   Plot and export all solve-fors in a filter or smoother run
#   (JSON input version)
#
#   Reads a GMAT filter or smoother JSON output file and plots each estimated
#   parameter beyond the six orbital states together with its 3-sigma uncertainty.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the JSON file.
#
#   usage: plot_solve_fors_json.py [-h] [--no_plot] [--no_csv] jsonfile
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

    state_names = data['CartesianStateNames']
    nstates = len(state_names)

    if nstates <= 6:
        print('No solve-for parameters beyond orbital states found in', jsonfile)
        sys.exit()

    param_names = state_names[6:]
    epochs = []
    param_vals = [[] for _ in param_names]
    param_sigs = [[] for _ in param_names]

    for u in updates:
        state = u['State']
        cov = u['Covariance']
        epochs.append(parse_epoch(u['Epoch']))
        for k, idx in enumerate(range(6, nstates)):
            param_vals[k].append(state[idx])
            param_sigs[k].append(np.sqrt(cov[idx][idx]))

    return epochs, runtype, [param_names, param_vals, param_sigs]


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    t, runtype, [param_names, param_vals, param_sigs] = data

    for i, param in enumerate(param_names):

        vals = np.array(param_vals[i])
        sigs = np.array(param_sigs[i])

        f = plt.figure()
        f.suptitle(runtype + ' ' + param + ' and 3-Sigma Uncertainty')

        plt.plot_date(t, vals, label=param,
            xdate=True, linestyle='solid', marker=None)

        plus_sigma  = vals + 3 * sigs
        minus_sigma = vals - 3 * sigs

        plt.plot_date(t, plus_sigma, label=param + ' 3-sigma uncertainty',
            xdate=True, linestyle='dashed', marker=None, color='green')
        plt.plot_date(t, minus_sigma, label=None,
            xdate=True, linestyle='dashed', marker=None, color='green')

        plt.fill_between(t, plus_sigma, minus_sigma,
            color='seagreen', alpha=0.1, label=None)

        plt.rc('grid', linestyle='dashed', color='gray', alpha=0.5)
        plt.grid()
        plt.legend(loc='lower center', ncol=3)
        plt.xlabel('Time UTC')

        f.autofmt_xdate()

        if save_plot:

            outfile = os.path.join(outdir, runtype + '_' + param + '.pdf')
            f.savefig(outfile, bbox_inches='tight')

        if save_csv:

            outfile = os.path.join(outdir, runtype + '_' + param + '.csv')
            with open(outfile, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Time', param, 'Sigma-' + param])
                outarray = np.array([t, vals, sigs], dtype=object)
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
