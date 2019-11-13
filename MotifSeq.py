import os
import sys
import gzip
import io
import traceback
import argparse
import math
import numpy as np
import scipy.stats as st
import h5py
import scrappy
from mlpy import dtw_subsequence
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
# rcParams['figure.figsize'] = [12.0, 12.0]
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
# import matplotlib.image as mpimg


'''

    James M. Ferguson (j.ferguson[at]garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    MotifSeq - finding signal motifs in raw nanopore signal

    --------------------------------------------------------------------------------------
    version 0.0 - initial
    ...
    version 1.3.0 - scrappie imported, cleaned up vis, python3, med-MAD scaling,
                    scrapped adapter and segments, change models,



    TODO:
        - Move methods of data processing into yield based functions
        - ensure all args removed from functions
        - make callable from other scripts
        - Take any signal format based on headers
        - MultiFast5File support
        - make readID a string with decode()

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
    VERSION = "1.3.0"

    parser = MyParser(
        description="MotifSeq - the Ctrl+f for signal. Signal-level local alignment of sequence motifs")
    group = parser.add_mutually_exclusive_group()
    mods = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--f5f",
                       help="File list of fast5 paths")
    group.add_argument("-p", "--f5_path",
                       help="Fast5 top dir")
    group.add_argument("-s", "--signal",
                       help="Extracted signal file from SquigglePull")
    parser.add_argument("-l", "--scale", default="medmad", choices=["zscale", "medmad"],
                       help="scaling/normalisation factor to use")
    mods.add_argument("-i", "--fasta_input",
                        help="fasta file to be converted to simulated signal by scrappy")
    parser.add_argument("--scrappie_model", default="squiggle_r94", choices=['squiggle_r94', 'squiggle_r94_rna', 'squiggle_r10'],
                        help="model to use with fasta_input for conversion")
    mods.add_argument("-m", "--model",
                        help="custom multiline .tsv of signal to search for - see docs - name{tab}60{tab}435...")
    parser.add_argument("-x", "--sig_extract", action="store_true",
                        help="Extract signal of match")
    parser.add_argument("--slope", type=float, default=2.90,
                        help="[Experimental] slope")
    parser.add_argument("--intercept", type=float, default=-9.6,
                        help="[Experimental] intercept")
    parser.add_argument("--std_const", type=float, default=0.08468,
                        help="[Experimental] standard deviation constant")
    parser.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    parser.add_argument("-scale_hi", "--scale_hi", type=int, default=1200,
                        help="Upper limit for signal outlier scaling")
    parser.add_argument("-scale_low", "--scale_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    parser.add_argument("-V", "--version", action="store_true",
                        help="Print version information")
    parser.add_argument("--verbose", action="store_true",
                        help="engage higher level of verbosity for troubleshooting")
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # print metadata
    if args.version:
        sys.stderr.write("SquiggleKit MotifSeq: {}\n".format(VERSION))
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode active - dumping info to stderr\n")
        sys.stderr.write("SquiggleKit MotifSeq: {}\n".format(VERSION))
        sys.stderr.write("args: {}\n".format(args))

    if args.scale == "zscale":
        import sklearn.preprocessing


    sys.stderr.write("\n\n**********************************************************\n")
    sys.stderr.write("*  z-score, p-value, probability, etc. are based on      *\n")
    sys.stderr.write("*     preliminary experimental modeling only             *\n")
    sys.stderr.write("*                Use at own risk                         *\n")
    sys.stderr.write("**********************************************************\n\n\n")

    squig = []

    if args.model:
        # model, m_order, L = read_synth_model(args.model)
        model, m_order, L = read_bait_model(args.model)
    if args.fasta_input:
        model, m_order, L = convert_fasta(args.fasta_input, args.scrappie_model)
    if args.sig_extract:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("fast5", "readID", "model", "start", "end", "length", "distance_score", "model_mean", "model_stdev", "Z-score", "p-value", "hit_Probability", "normalised_signal"))
    else:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("fast5", "readID", "model", "start", "end", "length", "distance_score", "model_mean", "model_stdev", "Z-score", "p-value", "hit_Probability"))

    if args.f5f:
        # file list of fast5 files.
        # fast5_name\tquality_score
        # not using the second column atm
        if args.f5f.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.f5f, 'rt') as s:
            for l in s:
                l = l.strip('\n')
                l = l.split('\t')[0]
                path = l
                l = l.split('/')
                fast5 = l[-1]
                sig, read_ID = process_fast5(path)
                if not sig:
                    sys.stderr.write("Failed to extract signal: {} {}\n".format(path, fast5))
                    continue
                sig = np.array(sig, dtype=int)
                sig = scale_outliers(sig, args)
                if args.scale == "zscale":
                    sig = sklearn.preprocessing.scale(sig,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)
                elif args.scale == "medmad":
                    arr = np.ma.array(sig).compressed()
                    med = np.median(arr)
                    mad = np.median(np.abs(arr - med))
                    scaled_mad = mad * 1.4826
                    mad_sig = []
                    for i in sig:
                        mad_sig.append((i - med) / scaled_mad)
                    sig = np.array(mad_sig)
                else:
                    sys.stderr.write("unknown scale parameter: {}\n".format(args.scale))
                    sys.exit()

                # Do the search
                if args.model:
                    get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)
                if args.fasta_input:
                    get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)


    elif args.f5_path:
        # process fast5 files given top level path
        for dirpath, dirnames, files in os.walk(args.f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)

                    # extract data from file
                    sig, read_ID = process_fast5(fast5_file)
                    if not sig:
                        sys.stderr.write("main():data not extracted. Moving to next file - {}\n".format(fast5_file))
                        continue
                    sig = np.array(sig, dtype=int)
                    sig = scale_outliers(sig, args)
                    if args.scale == "zscale":
                        sig = sklearn.preprocessing.scale(sig,
                                                      axis=0,
                                                      with_mean=True,
                                                      with_std=True,
                                                      copy=True)
                    elif args.scale == "medmad":
                        arr = np.ma.array(sig).compressed()
                        med = np.median(arr)
                        mad = np.median(np.abs(arr - med))
                        scaled_mad = mad * 1.4826
                        mad_sig = []
                        for i in sig:
                            mad_sig.append((i - med) / scaled_mad)
                        sig = np.array(mad_sig)
                    else:
                        sys.stderr.write("unknown scale parameter: {}\n".format(args.scale))
                        sys.exit()

                    # Do the search
                    if args.model:
                        get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)
                    if args.fasta_input:
                        get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)


    elif args.signal:
        # signal file, gzipped, from squigglepull
        # Header False for now, soon to be fixed - make argument, default True
        head = False
        if args.signal.endswith('.gz'):
            f_read = dicSwitch('gz')
        else:
            f_read = dicSwitch('norm')
        with f_read(args.signal, 'rt') as s:
            for l in s:
                if head:
                    head = False
                    continue
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                read_ID = l[1]
                # modify the l[6:] to the column the data starts...little bit of variability here.
                sig = np.array([float(i) for i in l[8:]])
                if not sig.any():
                    sys.stderr.write("No Signal found - please check signal format\n")
                    continue
                sig = scale_outliers(sig, args)
                if args.scale == "zscale":
                    sig = sklearn.preprocessing.scale(sig,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)
                elif args.scale == "medmad":
                    arr = np.ma.array(sig).compressed()
                    med = np.median(arr)
                    mad = np.median(np.abs(arr - med))
                    scaled_mad = mad * 1.4826
                    mad_sig = []
                    for i in sig:
                        mad_sig.append((i - med) / scaled_mad)
                    sig = np.array(mad_sig)
                else:
                    sys.stderr.write("unknown scale parameter: {}\n".format(args.scale))
                    sys.exit()

                # Do the search
                if args.model:
                    get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)
                if args.fasta_input:
                    get_region_multi(args, sig, model, m_order, fast5, read_ID, args.slope, args.intercept, args.std_const, L)

    else:
        sys.stderr.write("Unknown file or path input")
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


def scale_outliers(squig, args):
    '''
    Remove outliers based on hi/low args.
    I was scaling at one point, but removing tends to be less problematic
    This can change the position co-ordinates a little
    '''
    k = (squig > args.scale_low) & (squig < args.scale_hi)
    return squig[k]


def process_fast5(path):
    '''
    open fast5 and extract raw signal
    '''
    # open fast5 file
    squig = []
    readID = ""
    try:
        hdf = h5py.File(path, 'r')
    except:
        traceback.print_exc()
        sys.stderr.write("process_fast5():fast5 file failed to open: {}\n".format(path))
        squig = []
        return squig, readID
    # extract raw signal
    try:
        c = list(hdf['Raw/Reads'].keys())
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            squig.append(int(col))
        readID = hdf['Raw/Reads/'][c[0]].attrs['read_id']
    except:
        traceback.print_exc()
        sys.stderr.write("process_fast5():failed to extract events or fastq from: {}\n".format(path))
        squig = []
    return squig, readID


def read_synth_model(filename):
    '''
    read squiggle data from scrappie, old and new version, ready for dtw
    '''
    dic = {}
    m_order = []
    L = 0
    L_list = []
    with open(filename, 'rt') as r:
        for l in r:
            l = l.strip('\n')
            if l[0] == '#':
                if L != 0:
                    L_list.append(L)
                L = 0
                name = l[1:]
                dic[name] = []
                m_order.append(name)
            elif l[:3] == "pos":
                continue
            else:
                L += 1
                l = l.split()
                dic[name] = dic[name] + [float(l[2])] * int(round(float(l[4])))
        L_list.append(L)
    return dic, m_order, L_list


def convert_fasta(filename, mod):
    '''
    use scrappie to convert fasta sequences to dic of signals
    '''
    dic = {}
    m_order = []
    L_list = [] # lengths same order as m_order
    with open(filename, 'r') as r:
        for l in r:
            l = l.strip('\n')
            if l:
                if l[0] == ">":
                    name = l[1:]
                    dic[name] = []
                    # retain order from fasta file
                    m_order.append(name)
                    continue
                else:
                    L_list.append(len(l))
                    signal = scrappy.sequence_to_squiggle(l, model=mod).data(as_numpy=True, sloika=False)
                    for i in signal:
                        dwell =  int(round(math.exp(-i[2])))
                        dic[name] = dic[name] + [i[0]] * dwell
    return dic, m_order, L_list


def read_bait_model(filename):
    '''
    read baited model signal file - name \t kmer length \t signal....
    '''
    dic = {}
    m_order = []
    L_list = [] # lengths same order as m_order
    if filename.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(filename, 'rt') as s:
        for l in s:
            l = l.strip('\n')
            l = l.split('\t')
            name = l[0]
            L = int(l[1])
            # modify the l[3:] to the column the data starts...little bit of variability here.
            sig = np.array([float(i) for i in l[3:]], dtype=float)
            dic[name] = sig
    return dic, m_order, L_list


def get_region_multi(args, sig, model, m_order, fast5, readID, m, b, std, L):
    '''
    Find any region - simple demonstration version for 1 model
    '''
    c = 0
    for name in m_order:
        dist, cost, path = dtw_subsequence(model[name], sig)
        start = path[1][0]
        end = path[1][-1]
        # y = mx + b
        mod_mean = (m * L[c]) + b
        mod_stdev =  mod_mean * std
        Z = (dist - mod_mean) / mod_stdev
        p_value = st.norm.cdf(Z)
        hit_P = (1-p_value) * 100
        if args.sig_extract:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fast5, readID, name, start, end, end - start, dist, mod_mean, mod_stdev, Z, p_value, hit_P, '\t'.join([str(i) for i in sig[start:end]])))
        else:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fast5, readID, name, start, end, end - start, dist, mod_mean, mod_stdev, Z, p_value, hit_P))
        c += 1
        if args.view:
            view_region(sig, start, end, cost, path, model, dist, readID)
    return


def view_region(sig, start, end, cost, path, model, dist, readID):
    '''
    Visualise model position in Signal
    '''
    name = list(model.keys())[0]
    fig = plt.figure(1)
    fig.suptitle("readID: {}".format(readID), fontsize=16)

    ax1 = fig.add_subplot(221)   #top left
    ax1.axvline(x=start, color='m')
    ax1.axvline(x=end, color='m')

    plt.plot(sig, color='grey')
    plt.ylabel("Current (raw)")

    # fig = plt.figure(2)
    # ax = fig.add_subplot(111)
    # plot1 = plt.imshow(cost.T, origin='lower', cmap=cm.hot, interpolation='nearest')
    # plot2 = plt.plot(path[0], path[1], 'w')
    # xlim = ax2.set_xlim((-0.5, cost.shape[0]-0.5))
    # ylim = ax2.set_ylim((-0.5, cost.shape[1]-0.5))
    # # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    # plt.colorbar(cax=cax)

    # fig = plt.figure(3)
    ax2 = fig.add_subplot(222)   #top right
    # ax = fig.add_subplot(111)
    plt.plot(sig[start:end], color='grey')
    plt.plot(model[name], color='blue')
    grey_patch = mpatches.Patch(color='grey', label='raw signal')
    blue_patch = mpatches.Patch(color='blue', label='model signal')
    ax2.legend(handles=[grey_patch, blue_patch])

    # fig = plt.figure(4)
    # ax3 = fig.add_subplot(213)   #bottom left
    # ax = fig.add_subplot(111)
    # plt.plot(sig[start:end], color='grey')

    # fig = plt.figure(5)
    # ax = fig.add_subplot(111)
    # name = model.keys()[0]
    # plt.plot(model[name], color='blue')
    #
    # fig = plt.figure(6)
    ax3 = fig.add_subplot(223)   #bottom left
    plt.plot(cost[-1,])
    M = np.mean(cost[-1,])
    S = np.std(cost[-1,])
    bot = M - S
    ax3.axhline(y=dist, color='r')
    ax3.axhline(y=M, color='g')
    ax3.axhline(y=bot, color='b')
    plt.ylabel("distance score")
    green_patch = mpatches.Patch(color='green', label='mean')
    red_patch = mpatches.Patch(color='red', label='distance')
    blue_patch = mpatches.Patch(color='blue', label='1 stdev')
    ax3.legend(handles=[green_patch, red_patch, blue_patch])
    # ax4 = fig.add_subplot(224)
    # img=mpimg.imread('cappy.png')
    # plt.imshow(img)


    plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
