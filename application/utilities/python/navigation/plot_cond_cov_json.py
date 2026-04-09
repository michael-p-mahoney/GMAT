
#
#   Plot and export filter and smoother covariance condition numbers
#   (JSON input version)
#
#   Generates plots of the log (base 10) of the filter, backward filter, and
#   smoother covariance condition number.  The backward filter JSON file is
#   optional; if omitted only the forward filter and smoother curves are shown.
#
#   By default a PDF copy of the plot image and a CSV-file of the plot data are
#   saved in the same directory as the filter JSON file.
#
#   usage: plot_cond_cov_json.py [-h] [--no_plot] [--no_csv]
#                                filter_jsonfile smoother_jsonfile
#                                [--backfilter backfilter_jsonfile]
#
#   filter_jsonfile    (required) full path to filter JSON output file
#   smoother_jsonfile  (required) full path to smoother JSON output file
#   --backfilter       (optional) full path to backward filter JSON output file
#   --no_plot          (optional) don't show plot or create plot PDF files
#   --no_csv           (optional) don't create the CSV data files
#   -h, --help         (optional) show help message
#

import json
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, os, sys


def parse_epoch(s):
    """Parse a GMAT UTC epoch string like '10 Jun 2010 00:00:00.000'."""
    return datetime.strptime(s, "%d %b %Y %H:%M:%S.%f")


def load_updates(jsonfile):
    """Return (epochs, condition_numbers) for all updates in a JSON file."""

    with open(jsonfile) as f:
        data = json.load(f)

    if 'FilterUpdates' in data:
        updates = data['FilterUpdates']
    elif 'SmootherUpdates' in data:
        updates = data['SmootherUpdates']
    else:
        print('Error: Cannot find FilterUpdates or SmootherUpdates in', jsonfile)
        sys.exit()

    epochs = []
    cond   = []

    for u in updates:
        cov = np.array(u['Covariance'])
        if cov.size == 0:
            continue
        epochs.append(parse_epoch(u['Epoch']))
        cond.append(np.linalg.cond(cov))

    return epochs, np.array(cond)


def load(filter_jsonfile, smoother_jsonfile, backfilter_jsonfile=None):

    fil_t, fil_cond = load_updates(filter_jsonfile)
    smo_t, smo_cond = load_updates(smoother_jsonfile)

    bck_t, bck_cond = None, None
    if backfilter_jsonfile:
        bck_t, bck_cond = load_updates(backfilter_jsonfile)

    return [fil_t, fil_cond], [smo_t, smo_cond], [bck_t, bck_cond]


def render(data, outdir, show_plot=True, save_plot=True, save_csv=True):

    [fil_t, fil_cond], [smo_t, smo_cond], [bck_t, bck_cond] = data

    f = plt.figure()
    f.suptitle('log10(cond(covariance))')

    plt.plot_date(fil_t, np.log10(fil_cond), label='Forward Filter',
        xdate=True, linestyle='solid', fmt='C0')

    if bck_t is not None:
        plt.plot_date(bck_t, np.log10(bck_cond), label='Backward Filter',
            xdate=True, linestyle='solid', fmt='C1')

    plt.plot_date(smo_t, np.log10(smo_cond), label='Smoother',
        xdate=True, linestyle='solid', fmt='C2')

    plt.legend(loc='lower center', ncol=3, mode='expand')
    plt.xlabel('Time UTC')
    plt.ylabel('Unitless')

    f.autofmt_xdate()

    if save_plot:

        outfile = os.path.join(outdir, 'covariance_condition.pdf')
        f.savefig(outfile, bbox_inches='tight')

    if save_csv:

        outfile = os.path.join(outdir, 'covariance_condition.csv')
        with open(outfile, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            if bck_t is not None:
                csv_writer.writerow(['Time', 'Forward-Filter', 'Backward-Filter', 'Smoother'])
                outarray = np.array([fil_t, fil_cond, np.flip(bck_cond[1:]), smo_cond],
                                    dtype=object)
            else:
                csv_writer.writerow(['Time', 'Forward-Filter', 'Smoother'])
                outarray = np.array([fil_t, fil_cond, smo_cond], dtype=object)

            for row in outarray.transpose():
                csv_writer.writerow(row)

    if show_plot:
        plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('filter_jsonfile',
        help='Name of GMAT filter JSON output file')
    parser.add_argument('smoother_jsonfile',
        help='Name of GMAT smoother JSON output file')
    parser.add_argument('--backfilter',
        help='Name of GMAT backward filter JSON output file (optional)',
        default=None)
    parser.add_argument('--no_plot',
        help="Don't show or save plot", action='store_false')
    parser.add_argument('--no_csv',
        help="Don't save data to CSV file", action='store_false')

    args = parser.parse_args()
    outdir, _ = os.path.split(args.filter_jsonfile)

    data = load(args.filter_jsonfile, args.smoother_jsonfile, args.backfilter)

    render(data, outdir, show_plot=args.no_plot,
        save_plot=args.no_plot, save_csv=args.no_csv)
