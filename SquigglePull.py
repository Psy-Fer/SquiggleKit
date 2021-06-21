import os
import sys
import argparse
import traceback
import numpy as np
import h5py
import time

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
        - Multi fast5 file support
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
        description="SquigglePull - extraction and (optional) conversion to pA of raw signal from Oxford Nanopore fast5 files")

    # arguments
    parser.add_argument("-p", "--path",
                        help="Top directory path of fast5 files")
    parser.add_argument("-t", "--type", action="store", default="auto", choices=["auto", "single", "multi"], help="Specify the type of files provided. Default is autodetection which enables a mix of single and multifast5 files.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    parser.add_argument("-r", "--raw_signal", action="store_true",
                        help="No conversion to pA, raw signal is extracted instead")
    parser.add_argument("-i", "--extra_info", action="store_true",
                        help="Print extra information used for signal conversion and in methylation calling - nanopolish/f5c")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode on. Starting timer.\n")
        start_time = time.time()

    if not os.path.isdir(args.path):
        sys.stderr.write("The provided path {} is not an existing directory.\n".format(args.path))
        sys.exit(1)

    # process fast5 files given top level path
    # Changed this section to work with new function
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)
                # extract data from file
                data, multi = extract_f5_all(fast5_file, args)
                if not data:
                    sys.stderr.write("main():data not extracted from {}. Moving to next file.".format(fast5_file))
                    continue
                 print data
                if not multi:
                    print_data(data, args, fast5)
                else:
                    for read in data:
                        print_data(data[read], args, fast5)

    if args.verbose:
        end_time = time.time() - start_time
        sys.stderr.write("Time taken: {}\n".format(end_time))

# Added this function by combining the separate extraction functions for single and multi, and the pull function
def extract_f5_all(filename, args):
    '''
    inputs:
        filepath/name
        args from command line
    does:
        open fast5 files, extract whole signal and read data and converts to pA by default
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

#    with h5py.File(filename, 'r') as hdf:
#        multi = False
#        if args.type == "auto":
#            multi = f5_check_multi(hdf)
#            if args.verbose:
#                if multi:
#                    sys.stderr.write("{} detected as a multi fast5 file\n".format(filename))
#                else:
#                    sys.stderr.write("{} detected as a single fast5 file\n".format(filename)) 
#        elif args.type == "multi":
#            multi = True

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
                
                # convert to pA                    
                if not(args.raw_signal):
                    f5_dic['raw'] = np.array(f5_dic['raw'], dtype=int)
                    f5_dic['raw'] = convert_to_pA_numpy(f5_dic['raw'], digitisation, range, offset)
                    f5_dic['raw'] = np.round(f5_dic['raw'], 2)

                # save the extra info for printing                
                if args.extra_info:
                    f5_dic['digitisation'] = digitisation
                    f5_dic['offset'] = offset
                    f5_dic['range'] = range
                    f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5_all():failed to extract raw signal or fastq from {}".format(filename))
                f5_dic = {}

        # multi fast5 files
        else:
            for read in reads:
                f5_dic[read] = {'raw': [], 'seq': '', 'readID': '', 
                                'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                                'sampling_rate': 0.0}

                # extract the data
                try:
                    for col in hdf[read]['Raw/Signal'][()]:
                        f5_dic[read]['raw'].append(int(col))
                    f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                    digitisation = hdf[read]['channel_id'].attrs['digitisation']
                    offset = hdf[read]['channel_id'].attrs['offset']
                    range = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))

                    # convert to pA
                    if not(args.raw_signal):
                        f5_dic[read]['raw'] = np.array(f5_dic[read]['raw'], dtype=int)
                        f5_dic[read]['raw'] = convert_to_pA_numpy(f5_dic[read]['raw'], digitisation, range, offset)
                        f5_dic[read]['raw'] = np.round(f5_dic[read]['raw'], 2)
                    
                    # save the extra info for printing                    
                    if args.extra_info:
                        f5_dic[read]['digitisation'] = digitisation
                        f5_dic[read]['offset'] = offset
                        f5_dic[read]['range'] = range
                        f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
                    
                except:
                    traceback.print_exc()
                    sys.stderr.write("extract_fast5_all():failed to read readID: {}".format(read))
    
    return f5_dic, multi

# new numpy version of convert function
def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

# new function created to reduce redundancy of code
def print_data(data, args, fast5):

    ar = map(str, data['raw'])

    if args.extra_info:
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(fast5, data['readID'],
                data['digitisation'], data['offset'], data['range'],
                data['sampling_rate'], '\t'.join(ar)))
    else:                        
        print('{}\t{}\t{}'.format(
                fast5, data['readID'], '\t'.join(ar)))

# check whether file provided is multi or single
# doesn't work as cannot get file version attribute
def f5_check_multi(hdf):
    ver = hdf.get("file_version", default=None)
    if ver is not None:
        version = ver.split('.')
        if len(version) != 2:
            sys.write.stderr("Could not auto detect the file type. Please supply the type as an argument command line and re-run.")
            sys.exit(1)
        return version[0] >= 1
    else:
        return False

if __name__ == '__main__':
    main()

