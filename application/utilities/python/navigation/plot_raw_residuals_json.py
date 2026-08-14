
#
#   Plot filter raw residuals
#   (JSON input version)
#
#   Reads a GMAT filter JSON output file and plots the pre-update measurement
#   residuals, distinguishing accepted from sigma-edited observations.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the JSON file.
#
#   usage: plot_raw_residuals_json.py [-h] [--no_plot] [--no_csv] jsonfile
#
#   jsonfile    (required) full path to filter JSON output file
#   --no_plot   (optional) don't show plot or create plot PDF files
#   --no_csv    (optional) don't create the CSV data files
#   -h, --help  (optional) show help message
#

import json
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import argparse, os, sys


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

    rows = []
    for u in updates:
        if 'Residual' not in u:
            continue
        res = u['Residual']
        edit = u.get('EditFlag') or 'N'
        rows.append({
            'Time':       parse_epoch(u['Epoch']),
            'Residual-X': res[0] if len(res) > 0 else None,
            'Residual-Y': res[1] if len(res) > 1 else None,
            'Residual-Z': res[2] if len(res) > 2 else None,
            'EditFlag':   edit,
        })

    residuals = pd.DataFrame(rows)
    return residuals


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    residuals = data

    IACCEPT = residuals['EditFlag'].isin(['N'])
    ISIG    = residuals['EditFlag'].isin(['SIG'])

    t_accept = residuals.loc[IACCEPT, 'Time']
    t_sig    = residuals.loc[ISIG,    'Time']

    f = plt.figure()
    f.suptitle('Filter Raw Residuals')

    plt.plot(t_accept, residuals.loc[IACCEPT, 'Residual-X'], 'o', label='Residual-X')
    plt.plot(t_accept, residuals.loc[IACCEPT, 'Residual-Y'], 'o', label='Residual-Y')
    plt.plot(t_accept, residuals.loc[IACCEPT, 'Residual-Z'], 'o', label='Residual-Z')

    plt.plot(t_sig, residuals.loc[ISIG, 'Residual-X'], 'x', label='Residual-X-Edited', color='red')
    plt.plot(t_sig, residuals.loc[ISIG, 'Residual-Y'], 'x', label='Residual-Y-Edited', color='red')
    plt.plot(t_sig, residuals.loc[ISIG, 'Residual-Z'], 'x', label='Residual-Z-Edited', color='red')

    plt.legend(loc='lower center', ncol=4, mode='expand')
    plt.xlabel('Time UTC')
    plt.ylabel('Km')

    f.autofmt_xdate()

    if save_plot:

        outfile = os.path.join(outdir, 'filter_raw_residuals.pdf')
        f.savefig(outfile, bbox_inches='tight')

    if save_csv:

        outfile = os.path.join(outdir, 'filter_raw_residuals.csv')
        residuals.to_csv(outfile, index=False,
            columns=['Time', 'Residual-X', 'Residual-Y', 'Residual-Z', 'EditFlag'])

    if show_plot:
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('jsonfile',
        help='Name of GMAT filter JSON output file')
    parser.add_argument('--no_plot',
        help="Don't show or save plot", action='store_false')
    parser.add_argument('--no_csv',
        help="Don't save data to CSV file", action='store_false')

    args = parser.parse_args()
    outdir, _ = os.path.split(args.jsonfile)

    data = load(args.jsonfile)

    render(data, outdir, show_plot=args.no_plot,
        save_plot=args.no_plot, save_csv=args.no_csv)
