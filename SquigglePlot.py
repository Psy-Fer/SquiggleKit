import os
import sys
import gzip
import io
import traceback
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['figure.figsize'] = [16.0, 12.0]
matplotlib.rcParams['figure.dpi'] = 80
# matplotlib.rcParams['savefig.dpi'] = 300
import numpy as np
import h5py
# import sklearn.preprocessing
# import pandas as pd
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2019

    SquigglePlot - plotting the raw signal data and doing some basic manipulations

    --------------------------------------------------------------------------------------
    version 0.0 - initial


    TODO:
        - change size of saved figure
        - yaml/config file for multiple plots
        - Default plot settings rolled out tookit wide
        - light smoothing

    -----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2018 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    '''
    Main function for executing logic based on file input types
    '''
    parser = MyParser(
        description="SquigglePlot - plotting the raw signal data")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--f5f",
                       help="File list of fast5 paths")
    group.add_argument("-p", "--f5_path",
                       help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                       help="Extracted signal file from SquigglePull")
    group.add_argument("-i", "--ind",
                       help="Individual fast5 file")
    parser.add_argument("--head", action="store_true",
                       help="Header present in signal or flat file")
    parser.add_argument("-n", "--Num",
                        help="Section of signal to look at - -n 2000 or -n 100,1500")
    parser.add_argument("--scale_hi", type=int, default=1200,
                        help="Upper limit for signal outlier scaling")
    parser.add_argument("--scale_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    # Arguments for now, but best way forward will probably be a config file
    parser.add_argument("--plot_colour", default='grey',
                        help="Colour of signal plot, takes any pyplot entry: k,r,b,g,red,blue,etc...")
    parser.add_argument("--save",
                        help="Save file readname_saveArg.pdf --save saveArg.pdf, use png, etc for other file types")
    parser.add_argument("--save_path",
                        help="Save filepath")
    parser.add_argument("--no_show", action="store_true",
                        help="Do not show plot (used for saving many)")
    parser.add_argument("--dpi", type=int, default=100,
                        help="Change DPI for publication figs, eg: --dpi 300")


    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)


    matplotlib.rcParams['savefig.dpi'] = args.dpi

    N = 0
    N1 = 0
    N2 = 0
    if args.Num:
        if ',' in args.Num:
            N1, N2 = args.Num.split(',')
            N1, N2 = int(N1), int(N2)
        else:
            N = int(args.Num)

    head = False
    if args.head:
        head = True


    if args.f5f:
        # file list of fast5 files.
        # fast5_name\tquality_score
        # not using the second column atm
        if args.f5f.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.f5f, 'rb') as sz:
            if args.f5f.endswith('.gz'):
                sz = io.BufferedReader(sz)
            for l in sz:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')[0]
                path = l
                l = l.split('/')
                fast5 = l[-1]
                sig = process_fast5(path)
                if not sig:
                    continue
                if N:
                    sig = sig[:N]
                elif N1 or N2:
                    sig = sig[N1:N2]
                sig = scale_outliers(sig, args)
                # output sections
                view_sig(args, sig, fast5)

    elif args.f5_path:
        # process fast5 files given top level path
        for dirpath, dirnames, files in os.walk(args.f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)

                    # extract data from file
                    sig = process_fast5(fast5_file)
                    if not sig:
                        print >> sys.stderr, "main():data not extracted. Moving to next file", fast5_file
                        continue
                    if N:
                        sig = sig[:N]
                    elif N1 or N2:
                        sig = sig[N1:N2]
                    sig = np.array(sig, dtype=int)
                    sig = scale_outliers(sig, args)
                    view_sig(args, sig, fast5)

    elif args.signal:
        # signal file, gzipped, from squigglepull
        # testing
        if args.signal.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.signal, 'rb') as sz:
            if args.signal.endswith('.gz'):
                sz = io.BufferedReader(sz)
            for l in sz:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                # modify the l[6:] to the column the data starts...little bit of variability here.
                sig = np.array([int(i) for i in l[4:]], dtype=int)
                if not sig.any():
                    print >> sys.stderr, "nope 1"
                    continue
                if N:
                    sig = sig[:N]
                elif N1 or N2:
                    sig = sig[N1:N2]
                sig = scale_outliers(sig, args)
                view_sig(args, sig, fast5)

    elif args.ind:
        # Do an OS detection here for windows (get from fast5_fetcher)
        fast5 = args.ind.split('/')[-1]
        # extract data from file
        sig = process_fast5(args.ind)
        if not sig:
            print >> sys.stderr, "main():data not extracted.", args.ind
            parser.print_help(sys.stderr)
            sys.exit(1)
        if N:
            sig = sig[:N]
        elif N1 or N2:
            sig = sig[N1:N2]
        sig = np.array(sig, dtype=int)
        sig = scale_outliers(sig, args)
        view_sig(args, sig, fast5)

    else:
        print >> sys.stderr, "Unknown file or path input"
        parser.print_help(sys.stderr)
        sys.exit(1)

    print >> sys.stderr, "Done"


def dicSwitch(i):
    '''
    A switch to handle file opening and reduce duplicated code
    '''
    open_method = {
        "gz": gzip.open,
        "norm": open
    }
    return open_method[i]

def scale_outliers(sig, args):
    ''' Scale outliers to within m stdevs of median '''
    k = (sig > args.scale_low) & (sig < args.scale_hi)
    return sig[k]


def process_fast5(path):
    '''
    open fast5 and extract raw signal
    '''
    # open fast5 file
    sig = []
    try:
        hdf = h5py.File(path, 'r')
    except:
        traceback.print_exc()
        print >> sys.stderr, 'process_fast5():fast5 file failed to open: {}'.format(path)
        sig = []
        return sig
    # extract raw signal
    try:
        #b = sorted([i for i in hdf['Analyses'].keys() if i[0] == 'B'])[-1]
        c = hdf['Raw/Reads'].keys()
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            sig.append(int(col))
    except:
        traceback.print_exc()
        print >> sys.stderr, 'process_fast5():failed to extract events or fastq from', path
        sig = []
    return sig

def view_sig(args, sig, name, path=None):
    '''
    View the squiggle
    '''
    fig = plt.figure(1)
    # fig.subplots_adjust(hspace=0.1, wspace=0.01)
    # ax = fig.add_subplot(111)
    # plt.tight_layout()
    plt.autoscale()
    plt.title("Raw signal for:   {}".format(name))
    plt.xlabel("")
    plt.ylabel("Current (pA)")


    plt.plot(sig, color=args.plot_colour)
    if args.save:
        filename = os.path.join(args.save_path, "{}_dpi_{}_{}".format(name, args.dpi, args.save))
        plt.savefig(filename)
    if not args.no_show:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
