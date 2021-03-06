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

    Pull squiggle data from fast5 files

    input:
        - path to fast5 files

    output:
        - tsv signal file

    TODO:
        - Dynamic columns and data types
        - Mult fast5 file support
        - paf, sam, fastq, or flat file support for filtering
        - multiprocessing
        - use # notation at start of file for static values, size reduction,
          including things like kits, flowcells, versions, etc, for comparisons.



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
    parser.add_argument("--multi", action="store_true",
                        help="multi_fast5 files")
    parser.add_argument("-t", "--target",
                        help="Target information as comma delimited string structured by format type - SOON TO BE DEPRICATED")
    parser.add_argument("-f", "--form", default="all", choices=["pos1", "all"],
                        help="Format of target information - SOON TO BE DEPRICATED")
    group.add_argument("-r", "--raw", action="store_true",
                       help="Extract raw signal")
    group.add_argument("-e", "--event", action="store_true",
                       help="Extract event signal - SOON TO BE DEPRICATED")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    parser.add_argument("-s", "--scale", default=None, choices=["zscale", "medmad"],
                       help="scaling/normalisation factor to use")
    parser.add_argument("-c", "--pA_convert", action="store_true",
                        help="Convert raw signal to pA, for comparisons")
    parser.add_argument("-m", "--meth", action="store_true",
                        help="Print extra information used in methylation calling - nanopolish/f5c")
    # parser.add_argument("-a", "--paf",
    #                     help="paf alignment file for nt approach - Benchmarking")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode on. Starting timer")
        start_time = time.time()

    # if args.scale == "zscale":
    #     from sklearn import preprocessing


    # process fast5 files given top level path
    # This should work for multi-fast5 too, push detect into extract_f5()
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)

                # extract data from file
                if args.multi:
                    data = read_multi_fast5(fast5_file)
                    for read in data:
                        if args.pA_convert:
                            # convert signal to pA
                            pA_sig = convert_to_pA(data[read])
                            if args.scale:
                                sig = np.array(pA_sig, dtype=float)
                                pA_sig = scale_data(args, sig)
                            ar = []
                            for i in pA_sig:
                                ar.append(str(i))
                            print('{}\t{}\t{}'.format(
                                    fast5, data[read]['readID'], '\t'.join(ar)))
                        else:
                            signal = data[read]['raw']
                            if args.scale:
                                sig = np.array(data[read]['raw'])
                                signal = scale_data(args, sig)
                            ar = []
                            for i in signal:
                                ar.append(str(i))
                            print('{}\t{}\t{}'.format(fast5, data[read]['readID'], '\t'.join(ar)))
                else:
                    data = extract_f5(fast5_file, args)
                    if not data:
                        sys.stderr.write("main():data not extracted. Moving to next file")
                        continue
                    region = pull_target(data, args)

                    if not region:
                        sys.stderr.write("main():Region not found. Moving to next file")
                        continue

                    if args.event:
                        ar = []
                        for i in region[2]:
                            ar.append(str(i))
                        print('{}\t{}\t{}\t{}\t{}'.format(
                            fast5, data['readID'], region[0], region[1], '\t'.join(ar)))
                    elif args.raw:
                        if args.pA_convert:
                            # convert signal to pA
                            pA_sig = convert_to_pA(data)
                            ar = []
                            for i in pA_sig:
                                ar.append(str(i))
                            print('{}\t{}\t{}\t{}\t{}'.format(
                                    fast5, data['readID'], region[0], region[1], '\t'.join(ar)))
                        elif args.meth:
                            ar = []
                            for i in region[2]:
                                ar.append(str(i))
                            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(fast5, data['readID'],
                                    data['digitisation'], data['offset'], data['range'],
                                    data['sampling_rate'], '\t'.join(ar)))
                        else:
                            ar = []
                            for i in region[2]:
                                ar.append(str(i))
                            print('{}\t{}\t{}\t{}\t{}'.format(
                                fast5, data['readID'], region[0], region[1], '\t'.join(ar)))
    if args.verbose:
        end_time = time.time() - start_time
        sys.stderr.write("Time taken: {}\n".format(end_time))


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

    f5_dic = {'raw': [], 'events': [], 'seq': '', 'readID': '',
              'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
              'sampling_rate': 0.0}

    # open fast5 file
    try:
        hdf = h5py.File(filename, 'r')
    except:
        traceback.print_exc()
        sys.stderr.write("extract_fast5():fast5 file failed to open: {}".format(filename))
        f5_dic = {}
        return f5_dic

    # extract event signal
    if args.event:
        try:
            b = sorted([i for i in list(hdf['Analyses'].keys()) if i[0] == 'B'])[-1]
            c = list(hdf['Raw/Reads'].keys())
            for col in hdf['Analyses'][b]['BaseCalled_template']['Events'][()]:
                f5_dic['events'].append(float(col[0]))

            fq = hdf['Analyses'][b]['BaseCalled_template']['Fastq'][()
                                                                    ].split('\n')
            f5_dic['seq'] = fq
            f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()

        except:
            traceback.print_exc()
            sys.stderr.write("extract_fast5():failed to extract events or fastq from {}".format(filename))
            f5_dic = {}

    # extract raw signal
    elif args.raw:
        if args.pA_convert or args.meth:
            try:
                c = list(hdf['Raw/Reads'].keys())
                for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                    f5_dic['raw'].append(int(col))


                f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
                f5_dic['digitisation'] = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
                f5_dic['offset'] = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
                f5_dic['range'] = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
                f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']

            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5():failed to extract events or fastq from {}".format(filename))
                f5_dic = {}

        else:
            try:
                c = list(hdf['Raw/Reads'].keys())
                for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                    f5_dic['raw'].append(int(col))

                f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()

            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5():failed to extract events or fastq from {}".format(filename))
                f5_dic = {}

    # signal flag not set
    else:
        sys.stderr.write("extract_fast5():Please choose 'raw' or 'events' for the signal flag.")

    return f5_dic

def read_multi_fast5(filename):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    with h5py.File(filename, 'r') as hdf:
        for read in list(hdf.keys()):
            f5_dic[read] = {'raw': [], 'events': [], 'seq': '', 'readID': '', 'digitisation': 0.0,
                            'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}


            try:
                for col in hdf[read]['Raw/Signal'][()]:
                    f5_dic[read]['raw'].append(int(col))

                f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                f5_dic[read]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
                f5_dic[read]['offset'] = hdf[read]['channel_id'].attrs['offset']
                f5_dic[read]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
                f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5():failed to read readID: {}".format(read))
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

    dicf5_dic = {'events': [], 'moves': [], 'seq': '', 'readID': ''}
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
            signal = scale_data(args, signal)

        # region.append(data['readID'])
        region.append(target)
        region.append(target_type)
        region.append(signal)

    elif target_type == 'all':
        if args.raw:
            signal = np.array(data['raw'])
        else:
            signal = np.array(data['events'])
        if args.scale:
            signal = scale_data(args, signal)
        target = str(len(signal))

        #region.append(data['readID'])
        region.append(target)
        region.append(target_type)
        region.append(signal)

    else:
        sys.stderr.write("pull_target():target_type not recognised: {}".format(target_type))
        return default

    if region:
        return region
    else:
        sys.stderr.write("pull_target():Something went wrong. Region not found")
        return default


def scale_data(args, sig):
    '''
    Scale shift and scale for comparisons
    '''
    # try:
        # scaled_data = sklearn.preprocessing.scale(data,
        #                                           axis=0,
        #                                           with_mean=True,
        #                                           with_std=True,
        #                                           copy=True)
    # except:
    #     traceback.print_exc()
    #     sys.stderr.write("scale_data():Something went wrong, failed to scale data")
    #     return 0
    try:
        if args.scale == "zscale":
            scaled_data = sklearn.preprocessing.scale(sig,
                                          axis=0,
                                          with_mean=True,
                                          with_std=True,
                                          copy=True)
        elif args.scale == "medmad":
            arr = np.ma.array(sig, dtype=float).compressed()
            med = np.median(arr)
            mad = np.median(np.abs(arr - med))
            scaled_mad = mad * 1.4826
            mad_sig = []
            for i in sig:
                mad_sig.append((i - med) / scaled_mad)
            scaled_data = np.array(mad_sig)
        else:
            sys.stderr.write("unknown scale parameter: {}\n".format(args.scale))
            sys.exit()
    except:
        traceback.print_exc()
        sys.stderr.write("scale_data():Something went wrong, failed to scale data")
        return 0

    return scaled_data


def convert_to_pA(d):
    '''
    convert raw signal data to pA using digitisation, offset, and range
    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        rawptr[j] = (rawptr[j] + offset) * raw_unit;
    }
    '''
    digitisation = d['digitisation']
    range = d['range']
    offset = d['offset']
    raw_unit = range / digitisation
    new_raw = []
    for i in d['raw']:
        j = (i + offset) * raw_unit
        new_raw.append("{0:.2f}".format(round(j,2)))
    return new_raw


if __name__ == '__main__':
    main()
