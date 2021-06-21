## SquigglePull.py
import os
import sys
import argparse
import traceback
import numpy as np
import h5py
import time
import sklearn.preprocessing

## SquigglePlot.py
#import os
#import sys
import gzip
import io
#import traceback
#import argparse
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
#import numpy as np
#import h5py


## segmenter.py
#import numpy as np
#import matplotlib
matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
#import sys
#import os
#import gzip
#import io
#import argparse
#import traceback
#import h5py
#import sklearn.preprocessing

import subprocess


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = MyParser(
        description="SquiggleKit - set of tools for processing and manipulating Oxford Nanopore raw signal data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subcommand = parser.add_subparsers(help='subcommand --help for help messages', dest="command")
    
    # main options
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    # sub-module for fast5_fetcher

    # sub-module for SquigglePull
    pull = subcommand.add_parser('pull', help='extraction and (optional) conversion to pA of signal', 
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # SquigglePull sub-module options
    pull.add_argument("-p", "--path",
                        help="Top directory path of fast5 files")
    pull.add_argument("-t", "--type", action="store", default="auto", choices=["auto", "single", "multi"], help="Specify the type of files provided. Default is autodetection which enables a mix of single and multifast5 files.")
    pull.add_argument("-r", "--raw_signal", action="store_true",
                        help="No conversion to pA, raw signal is extracted instead")
    pull.add_argument("-i", "--extra_info", action="store_true",
                        help="Print extra information used for signal conversion and in methylation calling - nanopolish/f5c")    
    
    # sub-module for SquigglePlot
    plot = subcommand.add_parser('plot', help='plot the signal data after (optional) conversion to pA',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # SquigglePlot sub-module options
    group = plot.add_mutually_exclusive_group()
    #   need to make it so that raw_signal flag cannot be applied with the signal file
    #   currently cannot support it unless output file contains digitisation, range and offest values (extra info mode)

    group.add_argument("-p", "--f5_path",
                        help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                        help="Extracted signal file from SquigglePull")
    group.add_argument("-i", "--ind", nargs="+", 
                        help="Individual fast5 file/s.")
    plot.add_argument("-r", "--readID",
                        help="Individual readID to extract from a multifast5 file")
    #plot.add_argument("--single",action="store_true",
    #                    help="single fast5 files")
    plot.add_argument("-t", "--type", action="store", default="auto", choices=["auto", "single", "multi"], help="Specify the type of files provided. Default is autodetection which enables a mix of single and multifast5 files.")    
    plot.add_argument("--head", action="store_true",
                        help="Header present in signal or flat file")
    plot.add_argument("--raw_signal",action="store_true",
                        help="Plot raw signal instead of converting to pA")
    plot.add_argument("-n", "--Num",
                        help="Section of signal to look at - -n 2000 or -n 100,1500")
    plot.add_argument("--lim_hi", type=int, default=1200,
                        help="Upper limit for signal outliers")
    plot.add_argument("--lim_low", type=int, default=0,
                        help="Lower limit for signal outliers")
    # Arguments for now, but best way forward will probably be a config file
    plot.add_argument("--plot_colour", default='grey',
                        help="Colour of signal plot, takes any pyplot entry: k,r,b,g,red,blue,etc...")
    plot.add_argument("--save",
                        help="Save file readname_saveArg.pdf --save saveArg.pdf, use png, etc for other file types")
    plot.add_argument("--save_path",
                        help="Save filepath")
    plot.add_argument("--no_show", action="store_true",
                        help="Do not show plot (used for saving many)")
    plot.add_argument("--dpi", type=int, default=100,
                        help="Change DPI for publication figs, eg: --dpi 300")
    
    # sub-module for segmenter
    segment = subcommand.add_parser('segmenter', help='find obvious structural regions in squiggle data',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # segmenter sub-module options
    group = segment.add_mutually_exclusive_group()
    group.add_argument("-f", "--f5f",
                       help="File list of fast5 paths")
    group.add_argument("-p", "--f5_path",
                       help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                       help="Extracted signal file from squigglePull")
    segment.add_argument("--single",action="store_true",
                        help="single fast5 files")
    segment.add_argument("-n", "--Num", type=int, default=0,
                        help="Section of signal to look at - default 0=all")
    segment.add_argument("-e", "--error", type=int, default=5,
                        help="Allowable error in segment algorithm")
    segment.add_argument("-c", "--corrector", type=int, default=50,
                        help="Window size for increasing total error correction - better long segment detection")
    segment.add_argument("-w", "--window", type=int, default=150,
                        help="Minimum segment window size to be detected")
    segment.add_argument("-d", "--seg_dist", type=int, default=50,
                        help="Maximum distance between 2 segments to be merged into 1")
    segment.add_argument("-t", "--std_scale", type=float, default=0.75,
                        help="Scale factor of STDev about median")
    segment.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    segment.add_argument("-g", "--gap", action="store_true",
                        help="Turn on gap distance for stall to polyTAil")
    segment.add_argument("-b", "--gap_dist", type=int, default=3000,
                        help="Maximum distance between stall and polyTAil segment - for 10X/dRNA")
    segment.add_argument("-k", "--stall", action="store_true",
                        help="Turn on stall detection - must be present")
    segment.add_argument("-u", "--test", action="store_true",
                        help="Run Tests")
    segment.add_argument("-l", "--stall_len", type=float, default=0.25,
                        help="Minimum percentage of minimum window segment for initial stall segment")
    segment.add_argument("-j", "--stall_start", type=int, default=300,
                        help="Maximum distance for start of stall segment to be detected")
    segment.add_argument("-lim_hi", "--lim_hi", type=int, default=900,
                        help="Upper limit for signal outlier scaling")
    segment.add_argument("-lim_low", "--lim_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    segment.add_argument("--raw_signal",action="store_true",
                        help="Plot raw signal instead of converting to pA")

    # sub-module for MotifSeq

    # sub-module for web application
    web = subcommand.add_parser('web', help='launch the web application version of SquiggleKit',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # collect args
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode on. Starting timer.\n")
        start_time = time.time()

    if args.command == "pull":
        if len(sys.argv) == 2:
            pull.print_help(sys.stderr)
            sys.exit(1)
        squigglePull(args)

    if args.command == "plot":
        if len(sys.argv) == 2:
            plot.print_help(sys.stderr)
            sys.exit(1)
        squigglePlot(args)        

    if args.command == "segmenter":
        if len(sys.argv) == 2:
            segment.print_help(sys.stderr)
            sys.exit(1)
        segmenter(args)

    if args.command == "web":
        os.system("python main.py")

    if args.verbose:
        end_time = time.time() - start_time
        sys.stderr.write("Time taken: {}\n".format(end_time))
    
    sys.stderr.write("Done\n")

def squigglePull(args):
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)
                # extract data from file
                data, multi = extract_f5_all(fast5_file, args)
                if not data:
                    sys.stderr.write("main():data not extracted from {}. Moving to next file.".format(fast5_file))
                    continue
                # print data
                if not multi:
                    print_data(data, args, fast5)
                else:
                    for read in data:
                        print_data(data[read], args, fast5)

def squigglePlot(args):

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
                    data, multi = extract_f5_all(fast5_file, args)
                    if not multi:
                        sig = data.get('raw')
                        if not sig:
                            sys.stderr.write("main():data not extracted. Moving to next file: {}\n".format(fast5_file))
                            continue
                        if N:
                            sig = sig[:N]
                        elif N1 or N2:
                            sig = sig[N1:N2]
                        sig = np.array(sig, dtype=float)
                        sig = scale_outliers(sig, args)
                        view_sig(args, sig, fast5)
                    else:
                        # extract data from file
                        for read in data:
                            sig = data[read].get('raw')
                            if not sig:
                                sys.stderr.write("main():data not extracted from read {}. Moving to next file: {}\n".format(read, fast5_file))
                                continue
                            if N:
                                sig = sig[:N]
                            elif N1 or N2:
                                sig = sig[N1:N2]
                            sig = np.array(sig, dtype=float)
                            sig = scale_outliers(sig, args)
                            view_sig(args, sig, read)

    elif args.ind:
        files = args.ind
        for fast5_file in files:
            # Do an OS detection here for windows (get from fast5_fetcher)
            fast5_file = args.ind.split('/')[-1]
            # extract data from file
            sig = None
            data, multi = extract_f5_all(fast5_file, args)
            if not multi:
                sig = data.get('raw')
                read = fast5_file
            else:
                if args.readID:
                # if readID is provided, only get data with matching readID
                # is this the same for all multifast5 files?
                    sig = data['read_'+args.readID].get('raw')
                    read = args.readID
                else:
                    for read in data:
                        sig = data[read].get('raw')
                        if not sig.any():
                            sys.stderr.write("main():data not extracted from: {}".format(args.ind))
                        if N:
                            sig = sig[:N]
                        elif N1 or N2:
                            sig = sig[N1:N2]
                        sig = np.array(sig, dtype=float)
                        sig = scale_outliers(sig, args)
                        view_sig(args, sig, read)
                    continue

        if not sig.any():
            sys.stderr.write("main():data not extracted from: {}".format(args.ind))
            parser.print_help(sys.stderr)
            sys.exit(1)

        if N:
            sig = sig[:N]
        elif N1 or N2:
            sig = sig[N1:N2]

        sig = np.array(sig, dtype=float)
        sig = scale_outliers(sig, args)
        view_sig(args, sig, read)

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
                view_sig(args, sig, readID, fast5)

def segmenter(args):

    if not args.Num:
        args.Num = -1

    if args.f5f:
        # file list of fast5 files.
        # fast5_name {Tab} quality_score
        # not using the second column atm
        args.single = True
        with open(args.f5f, 'r') as s:
            for l in s:
                l = l.strip('\n')
                l = l.split('\t')[0]
                path = l
                l = l.split('/')
                fast5 = l[-1]
                data, multi = extract_f5_all(path, args)
                sig = data.get('raw')
                if not sig.any():
                    sys.stderr.write("main():data not extracted. Moving to next file: {}".format(fast5))
                    continue
                # cut signal based on -n flag
                sig = sig[:args.Num]
                # This removes very large high and low peaks
                sig = scale_outliers(sig, args)
                # Do the segment detection
                segs = get_segs(sig, args)
                if not segs:
                    sys.stderr.write("no segments found: {}".format(fast5))
                    continue
                # run tests on segments based on user question
                if args.test:
                    segs = test_segs(segs, args)
                    if not segs:
                        continue
                # output sections
                out = []
                for i, j in segs:
                    out.append(str(i))
                    out.append(str(j))
                    output = ",".join(out)
                print("\t".join([fast5, output]))
                # visualise for parameter tuning
                if args.view:
                    view_segs(segs, sig, args)

    elif args.f5_path:
        # process fast5 files given top level path, recursive
        for dirpath, dirnames, files in os.walk(args.f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)
                    # extract data from file
                    data, multi = extract_f5_all(fast5_file, args)
                    if args.single:
                        # extract data from file
                        sig = data.get('raw')
                        if not sig.any():
                            sys.stderr.write("main():data not extracted. Moving to next file: {}".format(fast5))
                            continue
                        # cut signal based on -n flag
                        sig = sig[:args.Num]
                        sig = np.array(sig, dtype=float)
                        # This removes very large high and low peaks
                        sig = scale_outliers(sig, args)
                        # Do the segment detection
                        segs = get_segs(sig, args)
                        if not segs:
                            sys.stderr.write("no segments found: {}".format(fast5))
                            continue
                        # run tests on segments based on user question
                        if args.test:
                            segs = test_segs(segs, args)
                            if not segs:
                                continue
                        # output sections
                        out = []
                        for i, j in segs:
                            out.append(str(i))
                            out.append(str(j))
                            output = ",".join(out)
                        print("\t".join([fast5, output]))
                        # visualise for parameter tuning
                        if args.view:
                            view_segs(segs, sig, args)
                    else:
                        for read in data:
                            sig = data[read].get('raw')
                            sig = sig[:args.Num]
                            if not sig.any():
                                sys.stderr.write("main():data not extracted from read {}. Moving to next file: {}\n".format(read, fast5_file))
                                continue
                            sig = np.array(sig, dtype=float)
                            sig = scale_outliers(sig, args)
                            segs = get_segs(sig, args)
                            if not segs:
                                sys.stderr.write("no segments found: {}".format(fast5))
                                continue
                            # run tests on segments based on user question
                            if args.test:
                                segs = test_segs(segs, args)
                                if not segs:
                                    continue
                            # output sections
                            out = []
                            for i, j in segs:
                                out.append(str(i))
                                out.append(str(j))
                                output = ",".join(out)
                            print("\t".join([read, output]))
                            # visualise for parameter tuning
                            if args.view:
                                view_segs(segs, sig, args)

    elif args.signal:
        # signal file, gzipped, from squigglepull
        head = False
        if args.signal.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.signal, 'r') as s:
            if args.signal.endswith('.gz'):
                s = io.BufferedReader(s)
            for l in s:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                # modify the l[3:] to the column the data starts...little bit of variability here.
                # TODO: make this based on column header
                if "." in l[4]:
                    sig = np.array([float(i) for i in l[4:]], dtype=float)
                else:
                    sig = np.array([int(i) for i in l[4:]], dtype=int)
                #sig = np.array([int(i) for i in l[4:]], dtype=int)
                if not sig.any():
                    sys.stderr.write("No signal found in file: {} {}".format(args.signal, fast5))
                    continue
                # cut signal based on -n flag
                sig = sig[:args.Num]
                # This removes very large high and low peaks
                sig = scale_outliers(sig, args)
                # Do the segment detection
                segs = get_segs(sig, args)
                if not segs:
                    sys.stderr.write("no segments found: {}".format(fast5))
                    continue
                # run tests on segments based on user question
                if args.test:
                    segs = test_segs(segs, args)
                    if not segs:
                        sys.stderr.write("no segs for testing: {}".format(fast5))
                        continue
                # output sections
                out = []
                for i, j in segs:
                    out.append(str(i))
                    out.append(str(j))
                    output = ",".join(out)
                print("\t".join([fast5, output]))
                # visualise for parameter tuning
                if args.view:
                    view_segs(segs, sig, args)

def extract_f5_all(filename, args):
    '''
    inputs:
        filepath/name
        args from command line
    does:
        open fast5 files, extract whole signal and read data
        converts to pA by default
    Returns:
        dic for further processing/printing
    '''
    f5_dic = {}
    multi = False
    with h5py.File(filename, 'r') as hdf:
        if args.type == "auto":
            reads = list(hdf.keys())
            if 'read' not in reads[1]:
                if args.verbose:
                    sys.stderr.write("{} detected as a single fast5 file\n".format(filename)) 
                multi = False
            else:
                if args.verbose:
                    sys.stderr.write("{} detected as a multi fast5 file\n".format(filename))
                multi = True
        elif args.type == "multi":
            reads = list(hdf.keys())
            multi = True

        # single fast5 files
        if not multi:
            f5_dic = {'raw': [], 'seq': '', 'readID': '',
                    'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                    'sampling_rate': 0.0}
            # extract the data
            try:
                c = list(hdf['Raw/Reads'].keys())
                for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                    f5_dic['raw'].append(int(col))
                f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
                digitisation = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
                offset = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
                range = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
                
                # convert to pA and round                    
                if not(args.raw_signal):
                    f5_dic['raw'] = np.array(f5_dic['raw'], dtype=int)
                    f5_dic['raw'] = convert_to_pA_numpy(f5_dic['raw'], digitisation, range, offset)
                    f5_dic['raw'] = np.round(f5_dic['raw'], 2)

                # save the extra info for printing                
                if args.command == "pull":
                    if args.extra_info:
                        f5_dic['digitisation'] = digitisation
                        f5_dic['offset'] = offset
                        f5_dic['range'] = range
                        f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5_all():failed to extract data from {}\n".format(filename))
                f5_dic = {}

        # multi fast5 files
        else:
            for read in list(hdf.keys()):
                f5_dic[read] = {'raw': [], 'seq': '', 'readID': '', 
                                'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                                'sampling_rate': 0.0}

                # extract the data
                try:
                    # if readID is provided, only extract and save that information
                    readID = hdf[read]['Raw'].attrs['read_id'].decode()
                    if args.command == 'plot':
                        if args.readID:
                            if readID != args.readID:
                                continue

                    f5_dic[read]['readID'] = readID
                    for col in hdf[read]['Raw/Signal'][()]:
                        f5_dic[read]['raw'].append(int(col))

                    digitisation = hdf[read]['channel_id'].attrs['digitisation']
                    offset = hdf[read]['channel_id'].attrs['offset']
                    range = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))

                    # convert to pA and round
                    if not(args.raw_signal):
                        f5_dic[read]['raw'] = np.array(f5_dic[read]['raw'], dtype=int)
                        f5_dic[read]['raw'] = convert_to_pA_numpy(f5_dic[read]['raw'], digitisation, range, offset)
                        f5_dic[read]['raw'] = np.round(f5_dic[read]['raw'], 2)
                    
                    # save the extra info for printing                    
                    if args.command == "pull":                    
                        if args.extra_info:
                            f5_dic[read]['digitisation'] = digitisation
                            f5_dic[read]['offset'] = offset
                            f5_dic[read]['range'] = range
                            f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
                    
                except:
                    traceback.print_exc()
                    sys.stderr.write("extract_fast5_all():failed to read readID: {}\n".format(read))

    return f5_dic, multi

#def check_multi(file, signal):
#    if signal:
        #check if signal file is a multi
        

def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

def scale_outliers(sig, args):
    ''' Scale outliers to within m stdevs of median '''
    ''' Remove outliers that don't fit within the specified bounds '''
    k = (sig > args.lim_low) & (sig < args.lim_hi)
    return sig[k]

def dicSwitch(i):
    '''
    A switch to handle file opening and reduce duplicated code
    '''
    open_method = {
        "gz": gzip.open,
        "norm": open
    }
    return open_method[i]

def print_data(data, args, fast5):

    ar = map(str, data['raw'])
    #for i in data['raw']:
    #    ar.append(str(i))

    if args.extra_info:
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(fast5, data['readID'],
                data['digitisation'], data['offset'], data['range'],
                data['sampling_rate'], '\t'.join(ar)))
    else:                        
        print('{}\t{}\t{}'.format(
                fast5, data['readID'], '\t'.join(ar)))

def view_sig(args, sig, name, fast5, path=None):
    '''
    View the squiggle
    '''
    fig = plt.figure(1)
    # fig.subplots_adjust(hspace=0.1, wspace=0.01)
    # ax = fig.add_subplot(111)
    # plt.tight_layout()
    plt.autoscale()
    plt.xlabel("")
    # print(sig.dtype)
    # print(sig.dtype == float64)
    
    if args.signal:
        if sig.dtype == float:
            raw = False    
        elif sig.dtype == int:
            raw = True
    else:
        raw = args.raw_signal
    if raw:
        plt.title("Raw signal for:   {}\nFile:   {}".format(name, fast5))
        plt.ylabel("Current - Not scaled")
    else:
        plt.title("Signal for:   {}\nFile:   {}".format(name, fast5))
        plt.ylabel("Current (pA)")       


    plt.plot(sig, color=args.plot_colour)
    if args.save:
        filename = os.path.join(args.save_path, "{}_dpi_{}_{}".format(name, args.dpi, args.save))
        plt.savefig(filename)
    if not args.no_show:
        plt.show()
    plt.clf()

def get_segs(sig, args):
    '''
    Get segments from signal
    This works by running through the signal and finding regions that are above
    the bot and below the top parameters, with some error tollerance, for a
    minimum window of length.
    '''

    mn = sig.min()
    mx = sig.max()
    mean = np.mean(sig)
    median = np.median(sig)
    # use this with outlier rejection to fix stdev thresholds
    stdev = np.std(sig)
    top = median + (stdev * args.std_scale)
    bot = median - (stdev * args.std_scale)

    # parameter tuning visualisation
    # TODO: Put tuning plots here

    # this is the algo. Simple yet effective
    prev = False  # previous string
    err = 0       # total error
    prev_err = 0  # consecutive error
    c = 0         # counter
    w = args.corrector        # window to increase total error thresh
    seg_dist = args.seg_dist  # distance between 2 segs to be merged as one
    start = 0     # start pos
    end = 0       # end pos
    segs = []     # segments [(start, stop)]
    for i in range(len(sig)):
        a = sig[i]
        if a < top and a > bot: # If datapoint is within range
            if not prev:
                start = i
                prev = True
            c += 1 # increase counter
            w += 1 # increase window corrector count
            if prev_err:
                prev_err = 0
            if c >= args.window and c >= w and not c % w: # if current window longer than detect limit, and corrector, and is divisible by corrector
                err -= 1 # drop current error count by 1
        else:
            if prev and err < args.error:
                c += 1
                err += 1
                prev_err += 1
                if c >= args.window and c >= w and not c % w:
                    err -= 1
            elif prev and (c >= args.window or not segs and c >= args.window * args.stall_len):
                end = i - prev_err # go back to where error stretch began for accurate cutting
                prev = False
                if segs and start - segs[-1][1] < seg_dist: # if segs very close, merge them
                    segs[-1][1] = end
                else:
                    segs.append([start, end]) # save segment
                c = 0
                err = 0
                prev_err = 0
            elif prev:
                prev = False
                c = 0
                err = 0
                prev_err = 0
            else:
                continue

    if segs:
        return segs
    else:
        # print >> sys.stderr, "no segs found"
        return False

def view_segs(segs, sig, args):
    '''
    View the segments on the squiggle
    '''
    fig = plt.figure(1)
    #fig.subplots_adjust(hspace=0.1, wspace=0.01)
    ax = fig.add_subplot(111)

    # Show segment lines
    for i, j in segs:
        ax.axvspan(i, j, alpha=0.5, color='m')

    plt.plot(sig, color='k')
    plt.show()
    plt.clf()

if __name__ == '__main__':
    main()
