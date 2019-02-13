import os
import sys
import gzip
import io
import traceback
import argparse
import numpy as np
import h5py
import sklearn.preprocessing
from mlpy import dtw_subsequence
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.cm as cm

'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    SigTools - a nonopore raw signal toolkit

    --------------------------------------------------------------------------------------
    version 0.0 - initial



    TODO:
        - Move methods of data processing into yield based functions
        - ensure all args removed from functions
        - make callable from other scripts

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
    do the thing
    '''
    parser = MyParser(
        description="MotifSeq - the Ctrl+f for signal. Signal-level local alignment of sequence motifs.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--f5f",
                       help="File list of fast5 paths")
    group.add_argument("-p", "--f5_path",
                       help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                       help="Extracted signal file from SquigglePull")
    parser.add_argument("-a", "--adapt",
                        help="Adapter model file")
    parser.add_argument("-m", "--model",
                        help="Query model file")
    # group.add_argument("-g", "--get_baits", choices=["pick", "auto"],
    #                     help="Generate baits file")
    parser.add_argument("--segs",
                        help="segmenter file")
    parser.add_argument("-b", "--baits",
                        help="signal bait file")
    parser.add_argument("-t", "--dtw_thresh",
                        help="DTW distance threshold for match")
    parser.add_argument("-d", "--motif_dist",
                        help="max distance of adapter from start of signal")
    parser.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    parser.add_argument("-scale_hi", "--scale_hi", type=int, default=1200,
                        help="Upper limit for signal outlier scaling")
    parser.add_argument("-scale_low", "--scale_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    squig = []

    if args.segs:
        segments = get_segs(args.segs)
    if args.adapt:
        adapter = read_synth_model(args.adapt)
    if args.model:
        model = read_synth_model(args.model)
    if args.baits:
        baits = read_bait_model(args.baits)

    if args.f5f:
        # file list of fast5 files.
        # fast5_name\tquality_score
        # not using the second column atm
        if args.f5f.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.f5f, 'rb') as s:
            if args.f5f.endswith('.gz'):
                s = io.BufferedReader(s)
            for l in s:
                l = l.strip('\n')
                l = l.split('\t')[0]
                path = l
                l = l.split('/')
                fast5 = l[-1]
                sig = process_fast5(path)
                if not sig:
                    print >> sys.stderr, "Failed to extract signal", path, fast5
                    continue
                sig = np.array(sig, dtype=int)
                sig = scale_outliers(sig, args.scale_hi, args.scale_low)
                sig = sklearn.preprocessing.scale(sig,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)

                # Do the search
                if args.segs and args.adapt:
                    if fast5 in segments:
                        get_adapter_2(args, sig, adapter, segments[fast5])
                    else:
                        get_adapter(args, sig, adapter)
                elif args.adapt:
                    get_adapter(args, sig, adapter)
                if args.model:
                    get_region(args, sig, model)

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
                    sig = np.array(sig, dtype=int)
                    sig = scale_outliers(sig, args.scale_hi, args.scale_low)
                    sig = sklearn.preprocessing.scale(sig,
                                                      axis=0,
                                                      with_mean=True,
                                                      with_std=True,
                                                      copy=True)
                    # Do the search
                    if args.segs and args.adapt:
                        if fast5 in segments:
                            get_adapter_2(args, sig, adapter, segments[fast5])
                        else:
                            get_adapter(args, sig, adapter)
                    elif args.adapt:
                        get_adapter(args, sig, adapter)
                    if args.model:
                        get_region(args, sig, model)

    elif args.signal:
        # signal file, gzipped, from squigglepull
        head = True
        if args.signal.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.signal, 'rb') as s:
            if args.signal.endswith('.gz'):
                s = io.BufferedReader(s)
            for l in s:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                # modify the l[6:] to the column the data starts...little bit of variability here.
                sig = np.array([float(i) for i in l[8:]])
                if not sig.any():
                    print >> sys.stderr, "nope 1"
                    continue
                sig = scale_outliers(sig, args.scale_hi, args.scale_low)
                sig = sklearn.preprocessing.scale(sig,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)
                # Do the search

                if args.segs and args.adapt:
                    if fast5 in segments:
                        get_adapter_2(args, sig, adapter, segments[fast5])
                    else:
                        get_adapter(args, sig, adapter)
                elif args.adapt:
                    get_adapter(args, sig, adapter)
                if args.model:
                    get_region(args, sig, model)

    else:
        print >> sys.stderr, "Unknown file or path input"
        parser.print_help(sys.stderr)
        sys.exit(1)


def dicSwitch(i):
    '''
    A switch to handle file opening and reduce duplicated code
    '''
    open_method = {
        "gz": gzip.open,
        "norm": open
    }
    return open_method[i]


def scale_outliers(squig, hi, low):
    ''' Scale outliers to within m stdevs of median '''
    ret = []
    for i in squig:
        if i > hi or i < low:
            continue
        ret.append(i)
    return np.array(ret, dtype=int)


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
        print >> sys.stderr, 'process_fast5():fast5 file failed to open: {}'.format(path)
        squig = []
        return squig
    # extract raw signal
    try:
        b = sorted([i for i in hdf['Analyses'].keys() if i[0] == 'B'])[-1]
        c = hdf['Raw/Reads'].keys()
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            squig.append(int(col))
    except:
        traceback.print_exc()
        print >> sys.stderr, 'process_fast5():failed to extract events or fastq from', path
        squig = []
    return squig


def read_synth_model(filename):
    '''
    read squiggle data ready for dtw
    '''
    dic = {}
    first = False
    with open(filename, 'r') as r:
        for l in r:
            l = l.strip('\n')
            if l[0] == '#':
                name = l[1:]
                dic[name] = []
            else:
                l = l.split()
                dic[name] = dic[name] + [float(l[2])] * int(round(float(l[4])))
    # print >> sys.stderr, len(dic['adapter'])
    return dic


def read_bait_model(filename):
    '''
    read baited signal file
    '''
    dic = {}
    if filename.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(filename, 'rb') as s:
        if filename.endswith('.gz'):
            s = io.BufferedReader(s)
        for l in s:
            l = l.strip('\n')
            l = l.split('\t')
            name = l[0]
            # modify the l[6:] to the column the data starts...little bit of variability here.
            sig = np.array([float(i) for i in l[3:]], dtype=float)
            dic[name] = sig
    return dic


def get_baits(args, sig, start, end, adapter_pos):
    """
    logic to get baits from candidate read by adjusting bountaries, and locking
    in a section to be pushed to a file
    """
    # loop on key press to modify boundaries and Visualise

    # exit method

    # accepty mods and return output

    return


def get_segs(segfile):
    '''
    create segment dic
    '''
    dic = {}
    with open(segfile, 'r') as f2:
        for l in f2:
            l = l.strip('\n')
            l = l.split()
            segs = l[1].split(',')
            dic[l[0]] = [int(segs[0]), int(segs[1])]
    return dic


def get_adapter(args, sig, adapter):
    '''
    Find adapter in signal using model/bait
    Call segmenter for Stall?
    '''
    pos = []
    sig_search = sig[:800]
    # seg = segmenter_call() # sig[seg[1]:]
    dist, cost, path = dtw_subsequence(adapter['adapter'], sig_search)
    start = path[1][0]
    end = path[1][-1]
    print "Dist:", dist, "pos:",  start, ",", end, "Dist from Start", start, "Length:", end - start
    if args.view:
        view_adapter(sig, start, end)

    if test_adapter(start):
        return True
    else:
        return False

    # pos = [dist, start, end]
    # return pos


def get_adapter_2(args, sig, adapter, segs):
    '''
    Find adapter in signal using model/bait
    Call segmenter for Stall?
    '''
    pos = []
    sig_search = sig[segs[1]:segs[1] + 800]
    # seg = segmenter_call() # sig[seg[1]:]
    dist, cost, path = dtw_subsequence(adapter['adapter'], sig_search)
    start = path[1][0] + segs[1]
    end = path[1][-1] + segs[1]

    print "Dist:", dist, "pos:",  start, ",", end, "Dist from Stall", start - \
        segs[1], "Length:", end - start
    if args.view:
        view_adapter(sig, start, end, s=segs)

    if test_adapter(start - segs[1]):
        return True
    else:
        return False
    #
    # pos = [dist, start, end]
    # return pos


def view_adapter(sig, start, end, s=False):
    '''
    Visualise adapter position in Signal
    '''

    fig = plt.figure(1)
    #fig.subplots_adjust(hspace=0.1, wspace=0.01)
    ax = fig.add_subplot(111)

    ax.axvline(x=start, color='m')
    ax.axvline(x=end, color='m')
    if s:
        ax.axvline(x=s[0], color='b')
        ax.axvline(x=s[1], color='b')
        ax.axvline(x=s[1] + 800, color='r')
    else:
        ax.axvline(x=800, color='r')

    plt.plot(sig, color='grey')
    plt.show()
    plt.clf()


def test_adapter(adapter_dist):
    '''
    look at positions of adapter and test distance from stall seg
    '''
    dist = 200
    if adapter_dist <= dist:
        return True
    else:
        return False


def get_region(args, sig, model):
    '''
    Find any region - simple demonstration version for 1 model
    '''
    pos = []
    name = model.keys()[0]
    # seg = segmenter_call() # sig[seg[1]:]
    dist, cost, path = dtw_subsequence(model[name], sig)
    start = path[1][0]
    end = path[1][-1]
    print "Dist:", dist, "pos:",  start, ",", end, "Dist from Start", start, "Length:", end - start
    if args.view:
        # view_region(sig, start, end)
        view_region(sig, start, end, cost, path, model)

    return

# def view_region(sig, start, end, s=False):


def view_region(sig, start, end, cost, path, model, s=False):
    '''
    Visualise model position in Signal
    '''

    fig = plt.figure(1)
    #fig.subplots_adjust(hspace=0.1, wspace=0.01)
    ax = fig.add_subplot(111)

    ax.axvline(x=start, color='m')
    ax.axvline(x=end, color='m')
    if s:
        ax.axvline(x=s[0], color='b')
        ax.axvline(x=s[1], color='b')

    plt.plot(sig, color='grey')

    # fig = plt.figure(2)
    # ax = fig.add_subplot(111)
    # plot1 = plt.imshow(cost.T, origin='lower', cmap=cm.hot, interpolation='nearest')
    # plot2 = plt.plot(path[0], path[1], 'w')
    # xlim = ax.set_xlim((-0.5, cost.shape[0]-0.5))
    # ylim = ax.set_ylim((-0.5, cost.shape[1]-0.5))
    # # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    # # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    # # plt.colorbar(cax=cax)

    # fig = plt.figure(3)
    # ax = fig.add_subplot(111)
    # plt.plot(sig[start:end], color='grey')
    # name = model.keys()[0]
    # plt.plot(model[name], color='blue')

    # fig = plt.figure(4)
    # ax = fig.add_subplot(111)
    # plt.plot(sig[start:end], color='grey')
    #
    # fig = plt.figure(5)
    # ax = fig.add_subplot(111)
    # name = model.keys()[0]
    # plt.plot(model[name], color='blue')

    plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
