import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import gzip
import argparse
import traceback
import h5py
import sklearn.preprocessing

'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    Segmenter - used to identify homopolymer/stall regions in signal data.

    --------------------------------------------------------------------------------------
    version 0.0 - initial



    TODO:
        - turn into a class to import and use easily
        - make yaml file for tuning args
        - Add parameter tuning args and plots


    -----------------------------------------------------------------------------
'''

def print_err(*args):
    sys.stderr.write(' '.join(map(str,args)) + '\n')

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    '''
    do the thing
    '''
    parser = MyParser(
        description="segmenter - script to find obvious regions in squiggle data")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--f5f",
                       help="File list of fast5 paths")
    group.add_argument("-p", "--f5_path",
                       help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                       help="Extracted signal file from squigglePull")
    parser.add_argument("-n", "--Num", type=int, default=0,
                        help="Section of signal to look at - default 0=all")
    parser.add_argument("-e", "--error", type=int, default=5,
                        help="Allowable error in segment algorithm")
    parser.add_argument("-c", "--corrector", type=int, default=50,
                        help="Window size for increasing total error correction - better long segment detection")
    parser.add_argument("-w", "--window", type=int, default=200,
                        help="Minimum segment window size to be detected")
    parser.add_argument("-d", "--seg_dist", type=int, default=50,
                        help="Maximum distance between 2 segments to be merged into 1")
    parser.add_argument("-t", "--std_scale", type=float, default=0.75,
                        help="Scale factor of STDev about median")
    parser.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    parser.add_argument("-g", "--gap", action="store_true",
                        help="Turn on gap distance for stall to polyTAil")
    parser.add_argument("-b", "--gap_dist", type=int, default=3000,
                        help="Maximum distance between stall and polyTAil segment - for 10X/dRNA")
    parser.add_argument("-k", "--stall", action="store_true",
                        help="Turn on stall detection - must be present")
    parser.add_argument("-u", "--test", action="store_true",
                        help="Run Tests")
    parser.add_argument("-l", "--stall_len", type=float, default=0.25,
                        help="Minimum percentage of minimum window segment for initial stall segment")
    parser.add_argument("-j", "--stall_start", type=int, default=300,
                        help="Maximum distance for start of stall segment to be detected")
    parser.add_argument("-scale_hi", "--scale_hi", type=int, default=900,
                        help="Upper limit for signal outlier scaling")
    parser.add_argument("-scale_low", "--scale_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not args.Num:
        args.Num = -1

    squig = []
    segs = []

    if args.f5f:
        # file list of fast5 files.
        # fast5_name {Tab} quality_score
        # not using the second column atm
        with open(args.f5f, 'r') as s:
            for l in s:
                l = l.strip('\n')
                l = l.split('\t')[0]
                path = l
                l = l.split('/')
                fast5 = l[-1]
                sig = process_fast5(path)
                if not sig:
                    continue
                sig = np.array(sig[:args.Num])
                sig = scale_outliers(sig, args)
                segs = get_segs(sig, args)
                if not segs:
                    continue
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

                if args.view:
                    view_segs(segs, sig, args)

    elif args.f5_path:
        # process fast5 files given top level path
        for dirpath, dirnames, files in os.walk(args.f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)

                    # extract data from file
                    sig = process_fast5(fast5_file)
                    if not sig:
                        print_err("main():data not extracted. Moving to next file", fast5_file)
                        continue
                    sig = sig[:args.Num]
                    sig = np.array(sig, dtype=int)
                    sig = scale_outliers(sig, args)
                    segs = get_segs(sig, args)
                    if not segs:
                        continue
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

                    if args.view:
                        view_segs(segs, sig, args)

    elif args.signal:
        # signal file, gzipped, from squigglepull
        head = True
        with gzip.open(args.signal, 'r') as s:
            for l in s:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                # modify the l[3:] to the column the data starts...little bit of variability here.
                sig = np.array([int(i) for i in l[3:]], dtype=int)
                if not sig.any():
                    print_err("nope 1")
                    continue
                sig = sig[:args.Num]
                sig = scale_outliers(sig, args)
                segs = get_segs(sig, args)
                if not segs:
                    print_err("nope 2")
                    continue
                if args.test:
                    segs = test_segs(segs, args)
                    if not segs:
                        print_err("nope 3")
                        continue
                # output sections
                out = []
                for i, j in segs:
                    out.append(str(i))
                    out.append(str(j))
                    output = ",".join(out)
                print("\t".join([fast5, output]))

                if args.view:
                    view_segs(segs, sig, args)

    else:
        print_err("Unknown file or path input")
        parser.print_help(sys.stderr)
        sys.exit(1)

    print_err("Done")


def scale_outliers(squig, args):
    ''' Scale outliers to within m stdevs of median '''
    k = (squig > args.scale_low) & (squig < args.scale_hi)
    return np.array(squig[k])


def process_fast5(path):
    '''
    open fast5 and extract raw signal
    '''
    # open fast5 file
    squig = []
    try:
        hdf = h5py.File(path, 'r')
    except:
        traceback.print_exc()
        print_err('process_fast5():fast5 file failed to open: {}'.format(path))
        squig = []
        return squig
    # extract raw signal
    try:
        #b = sorted([i for i in hdf['Analyses'].keys() if i[0] == 'B'])[-1]
        c = list(hdf['Raw/Reads'].keys())
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            squig.append(int(col))
    except:
        traceback.print_exc()
        print_err('process_fast5():failed to extract events or fastq from', path)
        squig = []
    return squig


def get_segs(sig, args):
    '''
    Get segments from signal
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
        if a < top and a > bot:
            if not prev:
                start = i
                prev = True
            c += 1
            w += 1
            if prev_err:
                prev_err = 0
            if c >= w and not c % w:
                err -= 1
        else:
            if prev and err < args.error:
                c += 1
                err += 1
                prev_err += 1
            elif prev and (c >= args.window or not segs and c >= args.window * args.stall_len):
                end = i - prev_err
                prev = False
                if segs and start - segs[-1][1] < seg_dist:
                    segs[-1][1] = end
                else:
                    segs.append([start, end])
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
        print_err("nope")
        return False


def test_segs(segs, args):
    '''
    test the segs meet various conditions
    ADD TESTS HERE!!!

    '''
    try:
        # Check that the first segement is close to beginning for stall
        if args.stall:
            if segs[0][0] > args.stall_start:
                print_err("start seg too late!")
                return False
        # Check second segment distance for polyT
        if args.gap:
            if segs[1][0] > segs[0][1] + args.gap_dist:
                print_err("second seg too far!!")
                return False
    except:
        print_err("nope!")
        traceback.print_exc()

    return segs


def view_segs(segs, sig, args):
    '''
    View the segments on the squiggle
    '''
    fig = plt.figure(1)
    #fig.subplots_adjust(hspace=0.1, wspace=0.01)
    ax = fig.add_subplot(111)

    for i, j in segs:
        ax.axvline(x=i, color='m')
        ax.axvline(x=j, color='m')

    plt.plot(sig, color='k')
    plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
