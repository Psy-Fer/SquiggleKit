import os
import sys
import argparse
import traceback
import numpy as np
import h5py
import time
import sklearn.preprocessing

'''
    SquigglePull
    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2017

    Pull squiggle data from fast5 files given some targetting method

    input:
        - path to fast5 files

    output:
        - tsv signal file

    TODO:
        - get meth parameters, for Hasindu

    // get channel parameters
    const char* scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    digitisation =
        fast5_read_float_attribute(scaling_group, "digitisation");
    offset = fast5_read_float_attribute(scaling_group, "offset");
    range = fast5_read_float_attribute(scaling_group, "range");
   sample_rate =fast5_read_float_attribute(scaling_group, "sampling_rate");


    Testing:
        python SquigglePull.py -es -p test/R9_event_data/ -t 20,110 -f pos1 > data.tsv
        python SquigglePull.py -r -p test/R9_raw_data/ -f all > data.tsv

    Notes:
        should do some target-type validations before executing and exit.
    -----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2017 James Ferguson

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
    One function to rule them all, one function to find them, One function to bring them
    all and in the darkness call them.
    '''

    parser = MyParser(
        description="SquigglePull - extraction of raw/event signal from Oxford Nanopore fast5 files")

    # arguments
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-p", "--path",
                        help="Top directory path of fast5 files")
    parser.add_argument("-t", "--target",
                        help="Target information as comma delimited string structured by format type")
    parser.add_argument("-f", "--form", default="all", choices=["pos1", "all"],
                        help="Format of target information")
    group.add_argument("-r", "--raw", action="store_true",
                       help="Target raw signal")
    group.add_argument("-e", "--event", action="store_true",
                       help="Target event signal")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    parser.add_argument("-s", "--scale", action="store_true",
                        help="Scale signal output for comparison")
    # parser.add_argument("-a", "--paf",
    #                     help="paf alignment file for nt approach - Benchmarking")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        print >> sys.stderr, "Verbose mode on. Starting timer"
        start_time = time.time()


    # process fast5 files given top level path
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)

                # extract data from file
                data = extract_f5(fast5_file, args)
                if not data:
                    print >> sys.stderr, "main():data not extracted. Moving to next file"
                    continue

                region = pull_target(data, args)

                if not region:
                    print >> sys.stderr, "main():Region not found. Moving to next file"
                    continue

                if args.event:
                    ar = []
                    for i in region[3]:
                        ar.append(str(i))
                    print '{}\t{}\t{}\t{}\t{}'.format(
                        fast5, region[0], region[1], region[2], '\t'.join(ar))
                elif args.raw:
                    ar = []
                    for i in region[2]:
                        ar.append(str(i))
                    print '{}\t{}\t{}\t{}'.format(
                    fast5, region[0], region[1], '\t'.join(ar))
    if args.verbose:
        end_time = time.time() - start_time
        print >> sys.stderr, "Time taken:", end_time


def extract_f5(filename, args, batch=False):
    '''
    inputs:
        filepath/name
        optional:
            Raw vs Events
            batch flags
    does:
        open fast5 files, extract whole signal and read data
    Returns:
        dic for further processing

    2 methods:
        - Open one at a time (slow) - single thread
        - Open batches at a time (fast) - paralellised


    Takes the latest basecalled events table.
    '''

    f5_dic = {'raw': [], 'events': [], 'seq': '', 'seq_id': ''}

    # open fast5 file
    try:
        hdf = h5py.File(filename, 'r')
    except:
        traceback.print_exc()
        print >> sys.stderr, 'extract_fast5():fast5 file failed to open: {}'.format(filename)
        f5_dic = {}
        return f5_dic

    # extract event signal
    if args.event:
        try:
            b = sorted([i for i in hdf['Analyses'].keys() if i[0] == 'B'])[-1]
            for col in hdf['Analyses'][b]['BaseCalled_template']['Events'][()]:
                f5_dic['events'].append(float(col[0]))

            fq = hdf['Analyses'][b]['BaseCalled_template']['Fastq'][()
                                                                    ].split('\n')
            f5_dic['seq'] = fq
            f5_dic['seq_id'] = fq[0].split(' ')[0].split('_')[0][1:]

        except:
            traceback.print_exc()
            print >> sys.stderr, 'extract_fast5():failed to extract events or fastq from', filename
            f5_dic = {}

    # extract raw signal
    elif args.raw:
        try:
            c = hdf['Raw/Reads'].keys()
            for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                f5_dic['raw'].append(int(col))

        except:
            traceback.print_exc()
            print >> sys.stderr, 'extract_fast5():failed to extract events or fastq from', filename
            f5_dic = {}

    # signal flag not set
    else:
        print >> sys.stderr, "extract_fast5():Please choose 'raw' or 'events' for the signal flag."

    return f5_dic


def pull_target(data, args, min_length=50, paf=None):
    '''
    Pull out target region from data.

    inputs:
        - data - dictionary containing reads
        - target - pos1: 20,110 - event/raw positions
        - target_type - pos1

    does:
        ...explain methods...

    Returns:
        - Regions of interest labelled by read_id/filename

    dicf5_dic = {'events': [], 'moves': [], 'seq': '', 'seq_id': ''}
    '''
    default = []
    region = []
    target_type = args.form
    if target_type not in ['all']:
        target = args.target.split(',')
        target = [int(i) for i in target]

    if target_type == 'pos1':
        # target: a,b
        if args.raw:
            signal = np.array(data['raw'][target[0]:target[1]])
        else:
            signal = np.array(data['events'][target[0]:target[1]])
        if args.scale:
            signal = scale_data(signal)

        region.append(data['seq_id'])
        region.append(target)
        region.append(target_type)
        region.append(signal)

    elif target_type == 'all':
        if args.raw:
            signal = np.array(data['raw'])
        else:
            signal = np.array(data['events'])
        if args.scale:
            signal = scale_data(signal)
        target = str(len(signal))

        #region.append(data['seq_id'])
        region.append(target)
        region.append(target_type)
        region.append(signal)

    else:
        print >> sys.stderr, "pull_target():target_type not recognised:", target_type
        return default

    if region:
        return region
    else:
        print >> sys.stderr, "pull_target():Something went wrong. Region not found"
        return default


def scale_data(data):
    '''
    Scale shift and scale for comparisons
    '''
    try:
        scaled_data = sklearn.preprocessing.scale(data,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)
    except:
        traceback.print_exc()
        print >> sys.stderr, "scale_data():Something went wrong, failed to scale data"
        return 0
    return scaled_data


if __name__ == '__main__':
    main()
