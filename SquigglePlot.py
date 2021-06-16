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
# matplotlib.rcParams['figure.figsize'] = [18.0, 12.0]
matplotlib.rcParams['figure.dpi'] = 80
# matplotlib.rcParams['savefig.dpi'] = 300
import numpy as np
import h5py

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
        - Update to handle pA scaled values with floats

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
        description="SquigglePlot - plotting the raw signal data after (optional) conversion to pA")
    group = parser.add_mutually_exclusive_group()
    #   need to make it so that raw_signal flag cannot be applied with the signal file
    #   currently cannot support it unless output file contains digitisation, range and offest values (extra info mode)

    group.add_argument("-p", "--f5_path",
                        help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                        help="Extracted signal file from SquigglePull. Currently not compatible with conversion")
    group.add_argument("-i", "--ind", nargs='+', 
                        help="Individual fast5 file/s")
    parser.add_argument("-r", "--readID",
                        help="Individual readID to extract from a multifast5 file")
    parser.add_argument("--single",action="store_true",
                        help="single fast5 files.")
    parser.add_argument("--head", action="store_true",
                        help="Header present in signal or flat file")
    parser.add_argument("--raw_signal",action="store_true",
                        help="Plot raw signal instead of converting to pA")
    parser.add_argument("-n", "--Num",
                        help="Section of signal to look at - -n 2000 or -n 100,1500")
    parser.add_argument("--lim_hi", type=int, default=1200,
                        help="Upper limit for signal outliers")
    parser.add_argument("--lim_low", type=int, default=0,
                        help="Lower limit for signal outliers")
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

    # check arguments are suitable and alert if not
    if args.raw_signal and args.signal:
        sys.stderr.write("Cannot convert the signal file. Will plot the values provided as is.\n")

    if args.save_path and not args.save:
        sys.stderr.write("Please provide a save suffix e.g. test.png using the --save argument.\n")
        sys.exit(1)

    if args.readID and args.single:
        sys.stderr.write("Run the program with the -i flag and supply only the single fast5 file with the read desired.\nIf you are unsure which fast5 file contains the desired read, use fast5fetcher to extract it.\n")
        sys.exit(1)

    if args.no_show and not args.save:
        sys.stderr.write("With the current settings, no output will be produced.\nEither remove the no_show flag, or select to save by providing a suffix to add to the created files using --save.\n")
        sys.exit(1)

    matplotlib.rcParams['savefig.dpi'] = args.dpi

    N = 0
    N1 = 0
    N2 = 0
    # Num gives section of the signal to look at
    if args.Num:
        if ',' in args.Num:
            N1, N2 = args.Num.split(',')
            N1, N2 = int(N1), int(N2)
        else:
            N = int(args.Num)

    head = False
    if args.head:
        head = True

    if args.f5_path:
        # process fast5 files given top level path
        for dirpath, dirnames, files in os.walk(args.f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)
                    #handle multifast5
                    if not args.single:
                        sigs = get_multi_fast5_signal(args, fast5_file)
                        for read in sigs:
                            sig = sigs[read]
                            if not sig.any():
                                sys.stderr.write("main():data not extracted from read {}. Moving to next file: {}\n".format(read, fast5_file))
                                continue
                            if N:
                                sig = sig[:N]
                            elif N1 or N2:
                                sig = sig[N1:N2]
                            sig = np.array(sig, dtype=float)
                            sig = scale_outliers(sig, args)
                            view_sig(args, sig, read, fast5_file)
                    else:
                        # extract data from file
                        sig = process_fast5(fast5_file, args)
                        if not sig.any():
                            sys.stderr.write("main():data not extracted. Moving to next file: {}".format(fast5_file))
                            continue
                        if N:
                            sig = sig[:N]
                        elif N1 or N2:
                            sig = sig[N1:N2]
                        sig = np.array(sig, dtype=float)
                        sig = scale_outliers(sig, args)
                        view_sig(args, sig, fast5, fast5_file)

    elif args.signal:
        # signal file, from squigglepull
        # testing
        with open(args.signal, 'rt') as sz:
            for l in sz:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                readID = l[1]
                if args.readID:
                    if args.readID != readID:
                        continue
                # modify the l[6:] to the column the data starts...little bit of variability here.
            # sig = np.array([int(i) for i in l[4:]], dtype=int)
                if "." in l[4]:
                    sig = np.array([float(i) for i in l[4:]], dtype=float)
                else:
                    sig = np.array([int(i) for i in l[4:]], dtype=int)
                if not sig.any():
                    sys.stderr.write("No signal found: {}".format(args.signal))
                    parser.print_help(sys.stderr)
                    sys.exit(1)
                if N:
                    sig = sig[:N]
                elif N1 or N2:
                    sig = sig[N1:N2]
                sig = scale_outliers(sig, args)
                if args.single:
                    view_sig(args, sig, fast5, fast5)
                else:
                    view_sig(args, sig, readID, fast5)
        
    elif args.ind:
        files = args.ind
        for fast5 in files:
        # Do an OS detection here for windows (get from fast5_fetcher)
            # extract data from file
            sig = None
            if args.single:
                sig = process_fast5(fast5, args)
                read = fast5.split('/')[-1]
                if not sig.any():
                    sys.stderr.write("main():data not extracted: {}".format(args.ind))
                    parser.print_help(sys.stderr)
                    sys.exit(1)
                if N:
                    sig = sig[:N]
                elif N1 or N2:
                    sig = sig[N1:N2]

                sig = np.array(sig, dtype=float)
                sig = scale_outliers(vimsig, args)
                view_sig(args, sig, read, fast5)

            else:
                sys.stderr.write("Looking at the file {}\n".format(fast5))
                sigs = get_multi_fast5_signal(args, fast5)
                if args.readID:
                # if readID is provided, only get data with matching readID
                    sig = sigs[args.readID]
                    read = args.readID
                    if not sig.any():
                        sys.stderr.write("main():data not extracted: {}".format(args.ind))
                        parser.print_help(sys.stderr)
                        sys.exit(1)
                    if N:
                        sig = sig[:N]
                    elif N1 or N2:
                        sig = sig[N1:N2]

                    sig = np.array(sig, dtype=float)
                    sig = scale_outliers(sig, args)
                    view_sig(args, sig, read, fast5)
                else:
                    for read in sigs:
                        sig = sigs[read]
                        if N:
                            sig = sig[:N]
                        elif N1 or N2:
                            sig = sig[N1:N2]
                        sig = np.array(sig, dtype=float)
                        sig = scale_outliers(sig, args)
                        view_sig(args, sig, read, fast5)            

    else:
        sys.stderr.write("Unknown file or path input")
        parser.print_help(sys.stderr)
        sys.exit(1)

    sys.stderr.write("Done\n")


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
    ''' Remove outliers that don't fit within the specified bounds '''
    k = (sig > args.lim_low) & (sig < args.lim_hi)
    return sig[k]

# Changed to convert the signal to pA
def get_multi_fast5_signal(args, read_filename):
    '''
    open multi fast5 files and extract information
    '''
    signals = {}
    f5_dic = read_multi_fast5(args, read_filename)
    for read in f5_dic:
        # get readID and signal
        readID = f5_dic[read]['readID']
        if args.readID:
            if readID != args.readID:
                continue
        signal = f5_dic[read]['signal']
        if not args.raw_signal:
            #convert to pA
            signal = np.array(signal, dtype=int)
            signal = convert_to_pA_numpy(signal, f5_dic[read]['digitisation'], f5_dic[read]['range'], f5_dic[read]['offset'])
            signal = np.round(signal, 2)
        signals[readID] = signal
    if args.readID and not signals:
        sys.stderr.write("Could not find data for read {} in file {}\n".format(args.readID, read_filename))
    # return signal/signals
    return signals

# Changed to only store information of certain read when readID is provided
# Also changed to extract information required for conversion
def read_multi_fast5(args, filename):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    with h5py.File(filename, 'r') as hdf:
        for read in list(hdf.keys()):
            f5_dic[read] = {'signal': [], 'readID': '', 'digitisation': 0.0,
                            'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}
            try:
                readID = hdf[read]['Raw'].attrs['read_id'].decode()
                #if readID is provided, only get and store data with matching readID
                if not args.readID is None:
                    if readID != args.readID:
                        continue
                f5_dic[read]['readID'] = readID
                f5_dic[read]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
                f5_dic[read]['offset'] = hdf[read]['channel_id'].attrs['offset']
                f5_dic[read]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
                f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']

                for col in hdf[read]['Raw/Signal'][()]:
                    f5_dic[read]['signal'].append(int(col))
            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5():failed to read readID: {}".format(read))
    return f5_dic

# Changed to extract information required for conversion and to convert the signal
def process_fast5(path, args):
    '''
    open fast5 and extract raw signal
    '''
    # open fast5 file
    sig = []
    try:
        hdf = h5py.File(path, 'r')
    except:
        traceback.print_exc()
        sys.stderr.write('process_fast5():fast5 file failed to open: {}'.format(path))
        sig = []
        return sig
    # extract raw signal
    try:
        #b = sorted([i for i in hdf['Analyses'].keys() if i[0] == 'B'])[-1]
        c = list(hdf['Raw/Reads'].keys())
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            sig.append(int(col))

        readID = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
        digitisation = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
        offset = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
        range = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
        if not args.raw_signal:
            #convert to pA
            sig = np.array(sig, dtype=int)
            sig = convert_to_pA_numpy(sig, digitisation, range, offset)
            sig = np.round(sig, 2)
            
    except:
        traceback.print_exc()
        sys.stderr.write('process_fast5():failed to extract events or fastq from: {}'.format(path))
        sig = []
    return sig

def view_sig(args, sig, name, file, path=None):
    '''
    View the squiggle
    '''
    fig = plt.figure(1)
    # fig.subplots_adjust(hspace=0.1, wspace=0.01)
    # ax = fig.add_subplot(111)
    # plt.tight_layout()
    plt.autoscale()
    plt.xlabel("")
    
    if args.signal:
        if sig.dtype == float:
            raw = False    
        elif sig.dtype == int:
            raw = True
    else:
        raw = args.raw_signal
    if raw:
        plt.title("Raw signal for:   {} \n File: {}".format(name, file))
        plt.ylabel("Current - Not scaled")
    else:
        plt.title("Signal for:   {} \nls File: {}".format(name, file))
        plt.ylabel("Current (pA)")       


    plt.plot(sig, color=args.plot_colour)
    if args.save:
        filename = os.path.join(args.save_path, "{}_dpi_{}_{}".format(name, args.dpi, args.save))
        plt.savefig(filename)
    if not args.no_show:
        plt.show()
    plt.clf()

# new conversion function (same as SquigglePull and segmenter)
def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

if __name__ == '__main__':
    main()
