
#
#   Plot filter scaled residuals
#   (JSON input version)
#
#   Reads a GMAT filter JSON output file and plots the scaled (normalised)
#   residuals.  Scaled residuals are available directly when the run used
#   DataFileStyle = 'Verbose'; otherwise they are not available and the script
#   will exit with an informative message.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the JSON file.
#
#   usage: plot_residual_ratios_json.py [-h] [--no_plot] [--no_csv] jsonfile
#
#   jsonfile    (required) full path to filter JSON output file
#   --no_plot   (optional) don't show plot or create plot PDF files
#   --no_csv    (optional) don't create the CSV data files
#   -h, --help  (optional) show help message
#

import json
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, os, sys


def parse_epoch(s):
    """Parse a GMAT UTC epoch string like '10 Jun 2010 00:00:00.000'."""
    return datetime.strptime(s, "%d %b %Y %H:%M:%S.%f")


def load(jsonfile):

    with open(jsonfile) as f:
        data = json.load(f)

    if 'FilterUpdates' in data:
        updates = data['FilterUpdates']
    elif 'SmootherUpdates' in data:
        updates = data['SmootherUpdates']
    else:
        print('Error: Cannot find FilterUpdates or SmootherUpdates in', jsonfile)
        sys.exit()

    # Collect only measurement updates that carry scaled residuals
    epochs = []
    scaled = []

    for u in updates:
        if 'ScaledResidual' not in u:
            continue
        res = u['ScaledResidual']
        epochs.append(parse_epoch(u['Epoch']))
        scaled.append(res)

    if not epochs:
        print('Error: No ScaledResidual data found in', jsonfile)
        print('Re-run GMAT with DataFileStyle = \'Verbose\' to include scaled residuals.')
        sys.exit()

    # Transpose: list-of-rows → per-component lists
    scaled = list(map(list, zip(*scaled)))   # scaled[i] = list of component-i values

    return epochs, scaled


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    t, scaled = data

    ncomp = len(scaled)
    labels = ['PosVec-X', 'PosVec-Y', 'PosVec-Z'] + \
             ['Component-' + str(i) for i in range(3, ncomp)]

    f = plt.figure()
    f.suptitle('Filter Scaled Residuals')

    for i in range(ncomp):
        lbl = labels[i] if i < len(labels) else 'Component-' + str(i)
        plt.plot(t, scaled[i], 'x', label=lbl)

    sigma = 3
    plus_sigma  = [+sigma] * len(t)
    minus_sigma = [-sigma] * len(t)

    plt.fill_between(t, plus_sigma, minus_sigma,
        color='seagreen', alpha=0.1, label='Accepted')

    plt.legend(loc='lower center', ncol=4, mode='expand')
    plt.xlabel('Time UTC')
    plt.ylabel('Unitless')

    f.autofmt_xdate()

    if save_plot:

        outfile = os.path.join(outdir, 'filter_scaled_residuals.pdf')
        f.savefig(outfile, bbox_inches='tight')

    if save_csv:

        outfile = os.path.join(outdir, 'filter_scaled_residuals.csv')
        with open(outfile, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            headers = ['Time'] + [labels[i] if i < len(labels) else 'Component-' + str(i)
                                   for i in range(ncomp)]
            csv_writer.writerow(headers)
            outarray = np.array([t] + scaled, dtype=object)
            for row in outarray.transpose():
                csv_writer.writerow(row)

    if show_plot:
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('jsonfile',
        help='Name of GMAT filter JSON output file (must use DataFileStyle = Verbose)')
    parser.add_argument('--no_plot',
        help="Don't show or save plot", action='store_false')
    parser.add_argument('--no_csv',
        help="Don't save data to CSV file", action='store_false')

    args = parser.parse_args()
    outdir, _ = os.path.split(args.jsonfile)

    data = load(args.jsonfile)

    render(data, outdir, show_plot=args.no_plot,
        save_plot=args.no_plot, save_csv=args.no_csv)
